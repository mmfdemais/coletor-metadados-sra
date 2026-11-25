import os
import sys
import pandas as pd

from coletor import normalizar, carregar_valores, resolver_valor

CHAVES_POSSIVEIS = ["Experiment", "Run", "BioSample", "Sample"]


def preparar(campo, valor, mapa_valores):
    if valor is None or (isinstance(valor, float) and pd.isna(valor)):
        return None
    if str(valor).strip() == "":
        return None
    final, _ = resolver_valor(campo, valor, mapa_valores)
    return final


def detectar_chave(df_manual, coluna_chave, df_ferramenta):
    valores = set(df_manual[coluna_chave].astype(str).str.strip())
    melhor, melhor_qtd = None, 0
    for coluna in CHAVES_POSSIVEIS:
        if coluna not in df_ferramenta.columns:
            continue
        qtd = len(valores & set(df_ferramenta[coluna].astype(str).str.strip()))
        if qtd > melhor_qtd:
            melhor, melhor_qtd = coluna, qtd
    return melhor, melhor_qtd


def comparar(caminho_manual, caminho_ferramenta, caminho_mapa,
             coluna_chave=None, caminho_saida="validacao.xlsx",
             aba_ferramenta="amostras", caminho_valores="valores.csv"):
    manual = pd.read_excel(caminho_manual).dropna(axis=1, how="all")
    ferramenta = pd.read_excel(caminho_ferramenta, sheet_name=aba_ferramenta)
    mapa = pd.read_csv(caminho_mapa)
    mapa_valores = carregar_valores(caminho_valores)

    if coluna_chave is None:
        coluna_chave = manual.columns[0]

    chave_ferramenta, casados = detectar_chave(manual, coluna_chave, ferramenta)
    if chave_ferramenta is None or casados == 0:
        raise ValueError("Nenhuma coluna de acesso em comum entre as planilhas.")

    manual = manual.copy()
    manual["__chave"] = manual[coluna_chave].astype(str).str.strip()
    manual = manual.drop_duplicates("__chave")

    ferramenta = ferramenta.copy()
    ferramenta["__chave"] = ferramenta[chave_ferramenta].astype(str).str.strip()
    ferramenta = ferramenta.drop_duplicates("__chave")

    juntas = manual.merge(ferramenta, on="__chave", how="inner",
                          suffixes=("_manual", "_ferramenta"))

    ausentes = manual[~manual["__chave"].isin(ferramenta["__chave"])]

    resumo = []
    divergencias = []

    for _, linha_mapa in mapa.iterrows():
        col_manual = str(linha_mapa["coluna_manual"]).strip()
        col_ferramenta = str(linha_mapa["campo_ferramenta"]).strip()

        if col_manual not in juntas.columns or col_ferramenta not in juntas.columns:
            resumo.append({"campo": col_ferramenta, "situacao": "coluna ausente",
                           "concordam": 0, "divergem": 0, "so_manual": 0,
                           "so_ferramenta": 0, "ambos_vazios": 0,
                           "cobertura_manual_pct": None, "acuracia_pct": None,
                           "indicios_novos": 0, "indicios_conferidos": 0,
                           "indicios_corretos": 0, "precisao_indicio_pct": None})
            continue

        concordam = divergem = so_manual = so_ferramenta = ambos_vazios = 0
        ganho_indicio = acerto_indicio = erro_indicio = 0
        coluna_indicio = f"{col_ferramenta}_indicio"
        tem_indicio = coluna_indicio in juntas.columns

        for _, linha in juntas.iterrows():
            bruto_m = preparar(col_ferramenta, linha[col_manual], mapa_valores)
            bruto_f = linha[col_ferramenta]
            tem_m = bruto_m is not None
            tem_f = pd.notna(bruto_f) and str(bruto_f).strip() != ""

            if not tem_f and tem_indicio:
                sugerido = preparar(col_ferramenta, linha[coluna_indicio],
                                    mapa_valores)
                if sugerido is not None:
                    ganho_indicio += 1
                    if tem_m:
                        alvo, achado = normalizar(bruto_m), normalizar(sugerido)
                        if achado and (achado in alvo or alvo in achado):
                            acerto_indicio += 1
                        else:
                            erro_indicio += 1

            if not tem_m and not tem_f:
                ambos_vazios += 1
            elif tem_m and not tem_f:
                so_manual += 1
                divergencias.append({"chave": linha["__chave"], "campo": col_ferramenta,
                                     "tipo": "so no manual",
                                     "manual": str(bruto_m).strip(), "ferramenta": "",
                                     "indicio": str(linha[coluna_indicio]).strip()
                                     if tem_indicio and pd.notna(linha[coluna_indicio]) else ""})
            elif tem_f and not tem_m:
                so_ferramenta += 1
                divergencias.append({"chave": linha["__chave"], "campo": col_ferramenta,
                                     "tipo": "so na ferramenta",
                                     "manual": "", "ferramenta": str(bruto_f).strip(),
                                     "indicio": ""})
            elif normalizar(bruto_m) == normalizar(bruto_f):
                concordam += 1
            else:
                divergem += 1
                divergencias.append({"chave": linha["__chave"], "campo": col_ferramenta,
                                     "tipo": "valor diferente",
                                     "manual": str(bruto_m).strip(),
                                     "ferramenta": str(bruto_f).strip(),
                                     "indicio": ""})

        preenchidos_manual = concordam + divergem + so_manual
        comparaveis = concordam + divergem
        conferiveis = acerto_indicio + erro_indicio

        resumo.append({
            "campo": col_ferramenta,
            "situacao": "comparado",
            "concordam": concordam,
            "divergem": divergem,
            "so_manual": so_manual,
            "so_ferramenta": so_ferramenta,
            "ambos_vazios": ambos_vazios,
            "cobertura_manual_pct": round(comparaveis / preenchidos_manual * 100, 1)
            if preenchidos_manual else None,
            "acuracia_pct": round(concordam / comparaveis * 100, 1)
            if comparaveis else None,
            "indicios_novos": ganho_indicio,
            "indicios_conferidos": conferiveis,
            "indicios_corretos": acerto_indicio,
            "precisao_indicio_pct": round(acerto_indicio / conferiveis * 100, 1)
            if conferiveis else None,
        })

    df_resumo = pd.DataFrame(resumo)
    df_divergencias = pd.DataFrame(divergencias)

    with pd.ExcelWriter(caminho_saida, engine="openpyxl") as writer:
        df_resumo.to_excel(writer, sheet_name="resumo", index=False)
        df_divergencias.to_excel(writer, sheet_name="divergencias", index=False)
        if not ausentes.empty:
            ausentes[[coluna_chave]].to_excel(
                writer, sheet_name="nao_recuperadas", index=False)

    return {
        "chave_manual": coluna_chave,
        "chave_ferramenta": chave_ferramenta,
        "manual_total": len(manual),
        "pareadas": len(juntas),
        "nao_recuperadas": len(ausentes),
        "resumo": df_resumo,
        "divergencias": df_divergencias,
    }


def main():
    if len(sys.argv) < 3:
        print("uso: python validar.py <planilha_manual.xlsx> <planilha_ferramenta.xlsx> "
              "[mapa_validacao.csv]")
        return

    caminho_mapa = sys.argv[3] if len(sys.argv) > 3 else "mapa_validacao.csv"
    if not os.path.exists(caminho_mapa):
        print(f"Mapa '{caminho_mapa}' nao encontrado.")
        return

    resultado = comparar(sys.argv[1], sys.argv[2], caminho_mapa)

    print(f"Chave: {resultado['chave_manual']} x {resultado['chave_ferramenta']}")
    print(f"Curadas manualmente: {resultado['manual_total']}")
    print(f"Pareadas: {resultado['pareadas']}")
    print(f"Nao recuperadas pela busca: {resultado['nao_recuperadas']}")
    print()
    print(resultado["resumo"].to_string(index=False))
    print()
    print("Detalhes em validacao.xlsx")


if __name__ == "__main__":
    main()
