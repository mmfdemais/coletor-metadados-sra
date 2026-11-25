from Bio import Entrez
import pandas as pd
import io
import os
import re
import time
import unicodedata
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

Entrez.email = "SEU_EMAIL_AQUI@exemplo.com"
Entrez.api_key = None

TERMO_BUSCA = '"B chromosome"[All Fields] AND fish[Organism]'
ARQUIVO_VOCABULARIO = "vocabulario.csv"
ARQUIVO_SAIDA = "amostras_ncbi.xlsx"
MANTER_VALOR_BRUTO = True
LIMITE_RESULTADOS = 10000


def normalizar(texto):
    if texto is None:
        return ""
    texto = unicodedata.normalize("NFKD", str(texto))
    texto = "".join(c for c in texto if not unicodedata.combining(c))
    texto = texto.lower().strip()
    texto = re.sub(r"[_\-/\.]+", " ", texto)
    texto = re.sub(r"[^a-z0-9 ]+", "", texto)
    texto = re.sub(r"\s+", " ", texto)
    return texto.strip()


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


def buscar_ids(termo, retmax=LIMITE_RESULTADOS):
    handle = Entrez.esearch(db="sra", term=termo, retmax=retmax)
    registro = Entrez.read(handle)
    handle.close()
    print(f"Encontrados {registro['Count']} resultados no SRA.")
    return registro["IdList"]


def buscar_runinfo(ids, tamanho_lote=200):
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
        print(f"  runinfo: {min(i + tamanho_lote, len(ids))}/{len(ids)}")
        time.sleep(0.4)
    return pd.concat(tabelas, ignore_index=True) if tabelas else pd.DataFrame()


def buscar_biosamples(biosample_ids, tamanho_lote=50):
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
            registros.append({
                "BioSample": amostra.get("accession"),
                "atributos": atributos,
            })
        print(f"  biosample: {min(i + tamanho_lote, len(biosample_ids))}/"
              f"{len(biosample_ids)}")
        time.sleep(0.4)
    return registros


def harmonizar(registros, mapa, ordem_campos):
    linhas = []
    nao_mapeados = Counter()
    exemplos = defaultdict(set)
    conflitos = []

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
                if MANTER_VALOR_BRUTO:
                    linha[f"{campo}__bruto"] = None
                continue
            valores = [v for _, v in entradas]
            distintos = list(dict.fromkeys(valores))
            linha[campo] = distintos[0]
            if MANTER_VALOR_BRUTO:
                linha[f"{campo}__bruto"] = " | ".join(
                    f"{n}={v}" for n, v in entradas)
            if len(distintos) > 1:
                conflitos.append({
                    "BioSample": registro["BioSample"],
                    "campo": campo,
                    "valores_encontrados": " | ".join(
                        f"{n}={v}" for n, v in entradas),
                })

        linhas.append(linha)

    df_nao_mapeados = pd.DataFrame([
        {
            "atributo_encontrado": nome,
            "ocorrencias": qtd,
            "exemplos_de_valor": " | ".join(sorted(exemplos[nome])),
            "campo_final_sugerido": "",
        }
        for nome, qtd in nao_mapeados.most_common()
    ])

    return pd.DataFrame(linhas), df_nao_mapeados, pd.DataFrame(conflitos)


def main():
    mapa, ordem_campos = carregar_vocabulario(ARQUIVO_VOCABULARIO)
    print(f"Vocabulario: {len(mapa)} sinonimos em "
          f"{len(ordem_campos)} campos finais.")

    ids = buscar_ids(TERMO_BUSCA)
    if not ids:
        print("Nenhum resultado. Revise o termo de busca.")
        return

    print("\nEtapa 1/3: metadados padronizados do SRA")
    df_runinfo = buscar_runinfo(ids)
    if "BioSample" not in df_runinfo.columns:
        print("Coluna BioSample ausente no runinfo.")
        return

    print("\nEtapa 2/3: atributos livres do BioSample")
    registros = buscar_biosamples(df_runinfo["BioSample"].dropna().tolist())

    print("\nEtapa 3/3: harmonizacao pelo vocabulario controlado")
    df_harmonizado, df_nao_mapeados, df_conflitos = harmonizar(
        registros, mapa, ordem_campos)

    colunas_sra = [c for c in ["Run", "Experiment", "Sample", "BioSample",
                               "BioProject", "ScientificName", "Platform",
                               "Model", "LibraryStrategy", "LibraryLayout",
                               "spots", "bases"]
                   if c in df_runinfo.columns]
    df_final = df_runinfo[colunas_sra].merge(
        df_harmonizado, on="BioSample", how="left")

    with pd.ExcelWriter(ARQUIVO_SAIDA, engine="openpyxl") as writer:
        df_final.to_excel(writer, sheet_name="amostras", index=False)
        df_nao_mapeados.to_excel(writer, sheet_name="nao_mapeados", index=False)
        df_conflitos.to_excel(writer, sheet_name="conflitos", index=False)

    total = len(df_final)
    print(f"\nConcluido: {total} linhas em '{ARQUIVO_SAIDA}'.")
    for campo in ordem_campos:
        if campo in df_final.columns:
            preenchidos = int(df_final[campo].notna().sum())
            pct = (preenchidos / total * 100) if total else 0
            print(f"  {campo}: {preenchidos}/{total} ({pct:.1f}%)")

    if not df_nao_mapeados.empty:
        print(f"\n{len(df_nao_mapeados)} atributos nao mapeados. "
              f"Veja a aba 'nao_mapeados', preencha 'campo_final_sugerido' "
              f"e adicione as linhas em '{ARQUIVO_VOCABULARIO}'.")


if __name__ == "__main__":
    main()
