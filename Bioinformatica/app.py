from flask import Flask, render_template, request, jsonify, send_file
import os
import threading
import traceback
import uuid
from datetime import datetime

import coletor

app = Flask(__name__)

PASTA_SAIDA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultados")
BASE = os.path.dirname(os.path.abspath(__file__))
ARQUIVO_VOCABULARIO = os.path.join(BASE, "vocabulario.csv")

CONFIGURACOES = {
    "vocabulario": (os.path.join(BASE, "vocabulario.csv"), "campo_final,sinonimo\n"),
    "valores": (os.path.join(BASE, "valores.csv"), "campo,valor_final,sinonimo\n"),
    "indicios": (os.path.join(BASE, "indicios.csv"), "campo,padrao\n"),
}

os.makedirs(PASTA_SAIDA, exist_ok=True)

tarefas = {}
trava = threading.Lock()


def registrar(tarefa_id, mensagem, tipo="info"):
    with trava:
        tarefa = tarefas.get(tarefa_id)
        if tarefa is None:
            return
        if tipo == "progresso" and tarefa["log"] and \
                tarefa["log"][-1]["tipo"] == "progresso":
            tarefa["log"][-1] = {"texto": mensagem, "tipo": tipo}
        else:
            tarefa["log"].append({"texto": mensagem, "tipo": tipo})


def rodar_coleta(tarefa_id, buscas, email, api_key, limite, manter_bruto,
                 buscar_texto):
    caminho = os.path.join(PASTA_SAIDA, f"amostras_{tarefa_id}.xlsx")
    try:
        resumo = coletor.coletar(
            buscas=buscas,
            email=email,
            caminho_vocabulario=CONFIGURACOES["vocabulario"][0],
            caminho_valores=CONFIGURACOES["valores"][0],
            caminho_indicios=CONFIGURACOES["indicios"][0],
            buscar_texto=buscar_texto,
            caminho_saida=caminho,
            api_key=api_key,
            limite=limite,
            manter_bruto=manter_bruto,
            log=lambda m, t="info": registrar(tarefa_id, m, t))
        with trava:
            tarefas[tarefa_id]["estado"] = "concluido"
            tarefas[tarefa_id]["resumo"] = resumo
            tarefas[tarefa_id]["caminho"] = caminho
        registrar(tarefa_id, f"{resumo['total']} amostras na planilha", "etapa")
    except Exception as erro:
        traceback.print_exc()
        with trava:
            tarefas[tarefa_id]["estado"] = "erro"
            tarefas[tarefa_id]["erro"] = str(erro)
        registrar(tarefa_id, str(erro), "erro")


@app.route("/")
def inicio():
    return render_template("index.html")


@app.route("/api/config/<nome>", methods=["GET", "POST"])
def configuracao(nome):
    if nome not in CONFIGURACOES:
        return jsonify({"erro": "Configuracao desconhecida."}), 404
    caminho, cabecalho = CONFIGURACOES[nome]

    if request.method == "POST":
        conteudo = request.json.get("conteudo", "")
        with open(caminho, "w", encoding="utf-8") as arquivo:
            arquivo.write(conteudo)
        return jsonify({"ok": True})

    if not os.path.exists(caminho):
        return jsonify({"conteudo": cabecalho})
    with open(caminho, encoding="utf-8") as arquivo:
        return jsonify({"conteudo": arquivo.read()})


@app.route("/api/coletar", methods=["POST"])
def iniciar_coleta():
    dados = request.json or {}
    buscas = [b.strip() for b in dados.get("buscas", []) if b.strip()]
    email = (dados.get("email") or "").strip()

    if not buscas:
        return jsonify({"erro": "Adicione ao menos uma busca."}), 400
    if "@" not in email:
        return jsonify({"erro": "O NCBI exige um e-mail valido."}), 400

    tarefa_id = uuid.uuid4().hex[:10]
    with trava:
        tarefas[tarefa_id] = {
            "estado": "rodando",
            "log": [],
            "resumo": None,
            "erro": None,
            "criado": datetime.now().isoformat(),
        }

    limite_bruto = dados.get("limite")
    try:
        limite = int(limite_bruto) if limite_bruto not in (None, "", 0, "0") else None
    except (TypeError, ValueError):
        limite = None

    thread = threading.Thread(
        target=rodar_coleta,
        args=(tarefa_id, buscas, email,
              (dados.get("api_key") or "").strip() or None,
              limite,
              bool(dados.get("manter_bruto")),
              bool(dados.get("buscar_texto"))),
        daemon=True)
    thread.start()

    return jsonify({"tarefa_id": tarefa_id})


@app.route("/api/status/<tarefa_id>")
def status(tarefa_id):
    with trava:
        tarefa = tarefas.get(tarefa_id)
        if tarefa is None:
            return jsonify({"erro": "Tarefa nao encontrada."}), 404
        return jsonify({
            "estado": tarefa["estado"],
            "log": tarefa["log"],
            "resumo": tarefa["resumo"],
            "erro": tarefa["erro"],
        })


@app.route("/api/download/<tarefa_id>")
def download(tarefa_id):
    with trava:
        tarefa = tarefas.get(tarefa_id)
    if not tarefa or tarefa.get("estado") != "concluido":
        return jsonify({"erro": "Planilha ainda nao esta pronta."}), 404
    return send_file(tarefa["caminho"], as_attachment=True,
                     download_name="amostras_ncbi.xlsx")


if __name__ == "__main__":
    print("Abra http://127.0.0.1:5000 no navegador")
    app.run(debug=False, port=5000)
