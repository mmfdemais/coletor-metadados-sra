from Bio import Entrez
import pandas as pd
import io
import os
import re
import time
import unicodedata
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

CAMPOS_SRA = ["Run", "Experiment", "Sample", "BioSample", "BioProject",
              "ScientificName", "Platform", "Model", "LibraryStrategy",
              "LibrarySource", "LibraryLayout", "LibraryName", "SampleName",
              "spots", "bases"]

JANELA_EVIDENCIA = 90

SENTINELA_VAZIO = "(vazio)"

NULOS = {
    "missing", "not collected", "not applicable", "not provided",
    "not determined", "unknown", "unspecified", "uncalculated",
    "na", "n a", "null", "restricted access", "not available",
}


def eh_nulo(valor):
    if valor is None:
        return True
    chave = normalizar(valor)
    return chave == "" or chave in NULOS


def normalizar(texto):
    if texto is None:
        return ""
    texto = unicodedata.normalize("NFKD", str(texto))
    texto = "".join(c for c in texto if not unicodedata.combining(c))
    texto = texto.lower().strip()
    texto = re.sub(r"[_\-/\.]+", " ", texto)
    texto = re.sub(r"[^a-z0-9 ]+", "", texto)
    return re.sub(r"\s+", " ", texto).strip()


def carregar_vocabulario(caminho):
    if not os.path.exists(caminho):
        raise FileNotFoundError(f"Arquivo '{caminho}' nao encontrado.")
    df = pd.read_csv(caminho)
    mapa = {}
    for _, linha in df.iterrows():
        campo = str(linha["campo_final"]).strip()
        sinonimo = normalizar(linha["sinonimo"])
        if sinonimo:
            mapa[sinonimo] = campo
    ordem = list(dict.fromkeys(str(c).strip() for c in df["campo_final"]))
    return mapa, ordem


def carregar_valores(caminho):
    if not caminho or not os.path.exists(caminho):
        return {}
    df = pd.read_csv(caminho)
    mapa = {}
    for _, linha in df.iterrows():
        campo = str(linha["campo"]).strip()
        chave = normalizar(linha["sinonimo"])
        if campo and chave:
            mapa[(campo, chave)] = str(linha["valor_final"]).strip()
    return mapa


def carregar_indicios(caminho):
    if not caminho or not os.path.exists(caminho):
        return []
    df = pd.read_csv(caminho)
    padroes = []
    for _, linha in df.iterrows():
        campo = str(linha["campo"]).strip()
        bruto = str(linha["padrao"]).strip()
        if not campo or not bruto:
            continue
        try:
            padroes.append((campo, re.compile(bruto, re.IGNORECASE)))
        except re.error:
            padroes.append((campo, re.compile(re.escape(bruto), re.IGNORECASE)))
    return padroes


def resolver_campo(nome_atributo, mapa):
    chave = normalizar(nome_atributo)
    if not chave:
        return None
    if chave in mapa:
        return mapa[chave]
    for sinonimo, campo in mapa.items():
        if len(sinonimo) >= 4 and sinonimo in chave:
            return campo
    return None


def resolver_valor(campo, valor, mapa_valores):
    if valor is None:
        return None, False
    chave = normalizar(valor)
    if chave and mapa_valores:
        final = mapa_valores.get((campo, chave))
        if final == SENTINELA_VAZIO:
            return None, False
        if final:
            return final, True
    if eh_nulo(valor):
        return None, False
    return valor, False


def _pausa(api_key):
    time.sleep(0.15 if api_key else 0.4)


def buscar_ids(termo, limite, api_key=None, log=None):
    handle = Entrez.esearch(db="sra", term=termo, retmax=0)
    registro = Entrez.read(handle)
    handle.close()
    total = int(registro["Count"])

    if total == 0:
        return [], 0

    alvo = total if not limite else min(limite, total)
    bloco = 5000
    ids = []

    while len(ids) < alvo:
        handle = Entrez.esearch(db="sra", term=termo,
                                retstart=len(ids),
                                retmax=min(bloco, alvo - len(ids)))
        parcial = Entrez.read(handle)
        handle.close()
        obtidos = parcial["IdList"]
        if not obtidos:
            break
        ids.extend(obtidos)
        if alvo > bloco and log:
            log(f"IDs {len(ids)}/{alvo}", "progresso")
        _pausa(api_key)

    return ids, total


def executar_buscas(buscas, limite, api_key, log):
    origem = defaultdict(list)
    ordem_ids = []
    for termo in buscas:
        try:
            ids, total = buscar_ids(termo, limite, api_key, log)
        except Exception as erro:
            log(f"Busca falhou: {termo} ({erro})", "erro")
            continue
        if limite and total > limite:
            log(f"{total} resultados para {termo}, coletando {len(ids)}", "aviso")
        else:
            log(f"{total} resultados para {termo}")
        for identificador in ids:
            if identificador not in origem:
                ordem_ids.append(identificador)
            origem[identificador].append(termo)
        _pausa(api_key)
    repetidos = sum(1 for v in origem.values() if len(v) > 1)
    if repetidos:
        log(f"{len(ordem_ids)} amostras unicas ({repetidos} em mais de uma busca)")
    else:
        log(f"{len(ordem_ids)} amostras unicas")
    return ordem_ids, origem


def estimar_minutos(quantidade, api_key, com_texto=False):
    pausa = 0.15 if api_key else 0.4
    lotes = (quantidade + 199) // 200 + (quantidade + 49) // 50
    if com_texto:
        lotes += max(1, quantidade // 40)
    return (lotes * (pausa + 1.2)) / 60


def buscar_runinfo(ids, api_key, log, tamanho_lote=200):
    tabelas = []
    for i in range(0, len(ids), tamanho_lote):
        lote = ids[i:i + tamanho_lote]
        handle = Entrez.efetch(db="sra", id=",".join(lote),
                               rettype="runinfo", retmode="text")
        conteudo = handle.read()
        handle.close()
        if isinstance(conteudo, bytes):
            conteudo = conteudo.decode("utf-8", errors="ignore")
        if conteudo.strip():
            tabelas.append(pd.read_csv(io.StringIO(conteudo)))
        log(f"SRA {min(i + tamanho_lote, len(ids))}/{len(ids)}", "progresso")
        _pausa(api_key)
    return pd.concat(tabelas, ignore_index=True) if tabelas else pd.DataFrame()


def buscar_biosamples(biosample_ids, api_key, log, tamanho_lote=50):
    biosample_ids = [b for b in dict.fromkeys(biosample_ids)
                     if isinstance(b, str) and b.strip()]
    registros = []
    for i in range(0, len(biosample_ids), tamanho_lote):
        lote = biosample_ids[i:i + tamanho_lote]
        handle = Entrez.efetch(db="biosample", id=",".join(lote), retmode="xml")
        xml_bruto = handle.read()
        handle.close()
        if isinstance(xml_bruto, bytes):
            xml_bruto = xml_bruto.decode("utf-8", errors="ignore")
        raiz = ET.fromstring(xml_bruto)
        for amostra in raiz.findall("BioSample"):
            atributos = {}
            bloco = amostra.find("Attributes")
            if bloco is not None:
                for attr in bloco.findall("Attribute"):
                    nome = (attr.get("attribute_name")
                            or attr.get("harmonized_name") or "")
                    valor = (attr.text or "").strip()
                    if nome and valor:
                        atributos[nome] = valor
            descricao = []
            for etiqueta in ("Title", "Paragraph", "Comment"):
                for no in amostra.iter(etiqueta):
                    if no.text and no.text.strip():
                        descricao.append(no.text.strip())
            registros.append({"BioSample": amostra.get("accession"),
                              "atributos": atributos,
                              "texto": " ".join(descricao)})
        log(f"BioSample {min(i + tamanho_lote, len(biosample_ids))}/"
            f"{len(biosample_ids)}", "progresso")
        _pausa(api_key)
    return registros


def buscar_titulos_sra(ids, api_key, log, tamanho_lote=100):
    titulos = {}
    for i in range(0, len(ids), tamanho_lote):
        lote = ids[i:i + tamanho_lote]
        try:
            handle = Entrez.efetch(db="sra", id=",".join(lote), retmode="xml")
            xml_bruto = handle.read()
            handle.close()
        except Exception as erro:
            log(f"Titulos do SRA falharam: {erro}", "erro")
            continue
        if isinstance(xml_bruto, bytes):
            xml_bruto = xml_bruto.decode("utf-8", errors="ignore")
        try:
            raiz = ET.fromstring(xml_bruto)
        except ET.ParseError:
            continue
        for pacote in raiz.iter("EXPERIMENT_PACKAGE"):
            acesso = None
            experimento = pacote.find("EXPERIMENT")
            if experimento is not None:
                acesso = experimento.get("accession")
            if not acesso:
                continue
            partes = []
            for etiqueta in ("TITLE", "LIBRARY_NAME", "DESIGN_DESCRIPTION",
                             "SCIENTIFIC_NAME", "VALUE", "TAG"):
                for no in pacote.iter(etiqueta):
                    if no.text and no.text.strip():
                        partes.append(no.text.strip())
            titulos[acesso] = " ".join(dict.fromkeys(partes))
        log(f"Titulos {min(i + tamanho_lote, len(ids))}/{len(ids)}", "progresso")
        _pausa(api_key)
    return titulos


def buscar_textos_bioproject(acessos, api_key, log, tamanho_lote=40):
    acessos = [a for a in dict.fromkeys(acessos)
               if isinstance(a, str) and a.strip()]
    textos = {}
    pubmed_por_projeto = defaultdict(list)

    for i in range(0, len(acessos), tamanho_lote):
        lote = acessos[i:i + tamanho_lote]
        try:
            handle = Entrez.esearch(db="bioproject", term=" OR ".join(lote),
                                    retmax=len(lote))
            encontrados = Entrez.read(handle)["IdList"]
            handle.close()
            _pausa(api_key)
            if not encontrados:
                continue
            handle = Entrez.efetch(db="bioproject", id=",".join(encontrados),
                                   retmode="xml")
            xml_bruto = handle.read()
            handle.close()
        except Exception as erro:
            log(f"BioProject falhou no lote {i // tamanho_lote + 1}: {erro}", "erro")
            continue

        if isinstance(xml_bruto, bytes):
            xml_bruto = xml_bruto.decode("utf-8", errors="ignore")
        try:
            raiz = ET.fromstring(xml_bruto)
        except ET.ParseError:
            continue

        for projeto in raiz.iter("DocumentSummary"):
            acesso = None
            for identificador in projeto.iter("ArchiveID"):
                acesso = identificador.get("accession")
                break
            if not acesso:
                continue
            partes = []
            for etiqueta in ("Title", "Name", "Description", "Objectives"):
                for no in projeto.iter(etiqueta):
                    if no.text and no.text.strip():
                        partes.append(no.text.strip())
            textos[acesso] = " ".join(partes)
            for publicacao in projeto.iter("Publication"):
                pmid = publicacao.get("id")
                if pmid and pmid.isdigit():
                    pubmed_por_projeto[acesso].append(pmid)

        log(f"BioProject {min(i + tamanho_lote, len(acessos))}/{len(acessos)}",
            "progresso")
        _pausa(api_key)

    return textos, pubmed_por_projeto


def buscar_resumos_pubmed(pubmed_por_projeto, api_key, log, tamanho_lote=100):
    todos = list(dict.fromkeys(
        p for lista in pubmed_por_projeto.values() for p in lista))
    if not todos:
        return {}

    resumos = {}
    for i in range(0, len(todos), tamanho_lote):
        lote = todos[i:i + tamanho_lote]
        try:
            handle = Entrez.efetch(db="pubmed", id=",".join(lote), retmode="xml")
            xml_bruto = handle.read()
            handle.close()
        except Exception as erro:
            log(f"PubMed falhou: {erro}", "erro")
            continue
        if isinstance(xml_bruto, bytes):
            xml_bruto = xml_bruto.decode("utf-8", errors="ignore")
        try:
            raiz = ET.fromstring(xml_bruto)
        except ET.ParseError:
            continue
        for artigo in raiz.iter("PubmedArticle"):
            pmid = None
            for no in artigo.iter("PMID"):
                pmid = no.text
                break
            partes = []
            for etiqueta in ("ArticleTitle", "AbstractText"):
                for no in artigo.iter(etiqueta):
                    if no.text and no.text.strip():
                        partes.append(no.text.strip())
            if pmid:
                resumos[pmid] = " ".join(partes)
        log(f"PubMed {min(i + tamanho_lote, len(todos))}/{len(todos)}",
            "progresso")
        _pausa(api_key)
    return resumos


def procurar_indicios(texto, padroes):
    achados = {}
    if not texto:
        return achados
    limpo = re.sub(r"[_/|]+", " ", str(texto))
    limpo = re.sub(r"\s+", " ", limpo).strip()
    for campo, expressao in padroes:
        if campo in achados:
            continue
        encontro = expressao.search(limpo)
        if not encontro:
            continue
        inicio = max(0, encontro.start() - JANELA_EVIDENCIA // 2)
        fim = min(len(limpo), encontro.end() + JANELA_EVIDENCIA // 2)
        trecho = limpo[inicio:fim].strip()
        if inicio > 0:
            trecho = "..." + trecho
        if fim < len(limpo):
            trecho = trecho + "..."
        achados[campo] = {"termo": encontro.group(0), "trecho": trecho}
    return achados


def harmonizar(registros, mapa, ordem_campos, manter_bruto, mapa_valores):
    linhas = []
    nao_mapeados = Counter()
    exemplos = defaultdict(set)
    conflitos = []
    valores_livres = Counter()
    normalizados = 0

    for registro in registros:
        linha = {"BioSample": registro["BioSample"]}
        acumulado = defaultdict(list)

        for nome, valor in registro["atributos"].items():
            campo = resolver_campo(nome, mapa)
            if campo is None:
                nao_mapeados[nome] += 1
                if len(exemplos[nome]) < 3:
                    exemplos[nome].add(valor)
                continue
            acumulado[campo].append((nome, valor))

        for campo in ordem_campos:
            entradas = acumulado.get(campo, [])
            if not entradas:
                linha[campo] = None
                if manter_bruto:
                    linha[f"{campo}__bruto"] = None
                continue
            distintos = list(dict.fromkeys(v for _, v in entradas))
            escolhido = distintos[0]
            final, foi_normalizado = resolver_valor(campo, escolhido, mapa_valores)
            if foi_normalizado:
                normalizados += 1
            else:
                valores_livres[(campo, escolhido)] += 1
            linha[campo] = final
            if manter_bruto:
                linha[f"{campo}__bruto"] = " | ".join(
                    f"{n}={v}" for n, v in entradas)
            if len(distintos) > 1:
                conflitos.append({
                    "BioSample": registro["BioSample"],
                    "campo": campo,
                    "valores_encontrados": " | ".join(
                        f"{n}={v}" for n, v in entradas)})

        linhas.append(linha)

    df_nao_mapeados = pd.DataFrame([
        {"atributo_encontrado": nome,
         "ocorrencias": qtd,
         "exemplos_de_valor": " | ".join(sorted(exemplos[nome])),
         "campo_final_sugerido": ""}
        for nome, qtd in nao_mapeados.most_common()])

    df_valores_livres = pd.DataFrame([
        {"campo": campo,
         "valor_encontrado": valor,
         "ocorrencias": qtd,
         "valor_final_sugerido": ""}
        for (campo, valor), qtd in valores_livres.most_common()])

    return (pd.DataFrame(linhas), df_nao_mapeados, pd.DataFrame(conflitos),
            df_valores_livres, normalizados)


def coletar(buscas, email, caminho_vocabulario, caminho_saida,
            api_key=None, limite=None, manter_bruto=True, log=None,
            caminho_valores=None, caminho_indicios=None, buscar_texto=False):
    if log is None:
        def log(mensagem, tipo="info"):
            print(mensagem)

    Entrez.email = email
    Entrez.api_key = api_key or None

    mapa, ordem_campos = carregar_vocabulario(caminho_vocabulario)
    mapa_valores = carregar_valores(caminho_valores)
    padroes = carregar_indicios(caminho_indicios) if buscar_texto else []

    log(f"Vocabulario: {len(mapa)} sinonimos em {len(ordem_campos)} campos")
    if mapa_valores:
        log(f"Valores: {len(mapa_valores)} equivalencias")
    if padroes:
        log(f"Indicios: {len(padroes)} padroes de texto")

    ids, origem = executar_buscas(buscas, limite, api_key, log)
    if not ids:
        raise ValueError("Nenhuma amostra encontrada. Revise os termos de busca.")

    minutos = estimar_minutos(len(ids), api_key, buscar_texto)
    if minutos >= 2:
        log(f"Tempo estimado: cerca de {round(minutos)} minutos", "aviso")

    log("Baixando metadados do SRA", "etapa")
    df_runinfo = buscar_runinfo(ids, api_key, log)
    if "BioSample" not in df_runinfo.columns:
        raise ValueError("O SRA nao retornou a coluna BioSample.")

    log("Baixando atributos do BioSample", "etapa")
    registros = buscar_biosamples(
        df_runinfo["BioSample"].dropna().tolist(), api_key, log)

    log("Harmonizando campos e valores", "etapa")
    (df_harmonizado, df_nao_mapeados, df_conflitos,
     df_valores_livres, normalizados) = harmonizar(
        registros, mapa, ordem_campos, manter_bruto, mapa_valores)

    colunas = [c for c in CAMPOS_SRA if c in df_runinfo.columns]
    df_final = df_runinfo[colunas].merge(
        df_harmonizado, on="BioSample", how="left")

    df_indicios = pd.DataFrame()
    if padroes and "BioProject" in df_final.columns:
        log("Procurando indicios no texto dos estudos", "etapa")
        titulos = buscar_titulos_sra(ids, api_key, log)
        acessos = df_final["BioProject"].dropna().unique().tolist()
        textos, pubmed_por_projeto = buscar_textos_bioproject(
            acessos, api_key, log)
        resumos = buscar_resumos_pubmed(pubmed_por_projeto, api_key, log)

        texto_amostra = {r["BioSample"]: r.get("texto", "") for r in registros}
        campos_indicio = list(dict.fromkeys(c for c, _ in padroes))
        linhas_indicio = []

        for _, linha in df_final.iterrows():
            projeto = linha.get("BioProject")
            fontes = []
            nomes = " ".join(str(linha.get(c)) for c in
                             ("LibraryName", "SampleName")
                             if pd.notna(linha.get(c)))
            titulo = titulos.get(str(linha.get("Experiment")), "")
            if titulo or nomes:
                fontes.append(("Titulo", f"{titulo} {nomes}".strip()))
            if texto_amostra.get(linha.get("BioSample")):
                fontes.append(("BioSample", texto_amostra[linha["BioSample"]]))
            if isinstance(projeto, str) and textos.get(projeto):
                fontes.append(("BioProject", textos[projeto]))
            for pmid in pubmed_por_projeto.get(projeto, []):
                if resumos.get(pmid):
                    fontes.append((f"PubMed {pmid}", resumos[pmid]))

            registro_indicio = {"BioSample": linha.get("BioSample")}
            achados_totais = {}
            for nome_fonte, texto in fontes:
                for campo, achado in procurar_indicios(texto, padroes).items():
                    if campo not in achados_totais:
                        achados_totais[campo] = (nome_fonte, achado)

            for campo in campos_indicio:
                if campo in achados_totais:
                    nome_fonte, achado = achados_totais[campo]
                    registro_indicio[f"{campo}_indicio"] = achado["termo"]
                    registro_indicio[f"{campo}_fonte"] = nome_fonte
                    registro_indicio[f"{campo}_evidencia"] = achado["trecho"]
                else:
                    registro_indicio[f"{campo}_indicio"] = None
                    registro_indicio[f"{campo}_fonte"] = None
                    registro_indicio[f"{campo}_evidencia"] = None
            linhas_indicio.append(registro_indicio)

        df_indicios = pd.DataFrame(linhas_indicio).drop_duplicates("BioSample")
        df_final = df_final.merge(df_indicios, on="BioSample", how="left")

    if len(buscas) > 1:
        chave = "Experiment" if "Experiment" in df_final.columns else "Run"
        mapa_origem = {k: " ; ".join(v) for k, v in origem.items()}
        df_final["busca_origem"] = df_final[chave].map(
            lambda x: mapa_origem.get(str(x), ""))

    with pd.ExcelWriter(caminho_saida, engine="openpyxl") as writer:
        df_final.to_excel(writer, sheet_name="amostras", index=False)
        df_nao_mapeados.to_excel(writer, sheet_name="campos_nao_mapeados",
                                 index=False)
        df_valores_livres.to_excel(writer, sheet_name="valores_nao_mapeados",
                                   index=False)
        df_conflitos.to_excel(writer, sheet_name="conflitos", index=False)

    total = len(df_final)
    cobertura = []
    for campo in ordem_campos:
        if campo in df_final.columns:
            preenchidos = int(df_final[campo].notna().sum())
            item = {"campo": campo,
                    "preenchidos": preenchidos,
                    "total": total,
                    "pct": round(preenchidos / total * 100, 1) if total else 0.0}
            coluna_indicio = f"{campo}_indicio"
            if coluna_indicio in df_final.columns:
                so_indicio = int((df_final[campo].isna() &
                                  df_final[coluna_indicio].notna()).sum())
                item["indicios"] = so_indicio
            cobertura.append(item)

    return {
        "total": total,
        "cobertura": cobertura,
        "nao_mapeados": df_nao_mapeados.to_dict(orient="records"),
        "valores_livres": df_valores_livres.head(30).to_dict(orient="records"),
        "conflitos": len(df_conflitos),
        "normalizados": normalizados,
        "arquivo": os.path.basename(caminho_saida),
    }
