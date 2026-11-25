# Coletor de metadados do NCBI SRA

Ferramenta para coleta automatizada e harmonização de metadados do NCBI Sequence Read Archive, voltada a organismos não-modelo, para os quais não existem ontologias consolidadas.

A harmonização é feita por um vocabulário controlado curado pelo próprio pesquisador, em arquivos de texto editáveis pela interface, sem alteração do código-fonte.

Desenvolvida como Trabalho de Conclusão de Curso em Biomedicina (FAIP, Marília – SP, 2026).

## O que faz

- Recupera registros do SRA e os atributos associados em BioSample, pelas E-utilities do NCBI.
- Reconcilia nomes de atributo divergentes (`pop`, `population`, `ecotype`) em um campo único.
- Padroniza valores (`male`, `Male`, `M`) e trata marcadores nulos (`missing`, `not collected`) como ausência real.
- Recupera informação não declarada nos atributos buscando padrões no título do experimento, na descrição do BioSample, no BioProject e no resumo do artigo vinculado no PubMed, registrando a fonte de cada indício.
- Entrega uma planilha `.xlsx` com uma linha por amostra.

## Requisitos

Python 3.10 ou superior.

```
pip install -r requirements.txt
```

## Uso

```
python app.py
```

Abra `http://127.0.0.1:5000` no navegador. Informe os termos de busca e um e-mail (exigido pelo NCBI). A chave de API é opcional e eleva o limite de requisições.

Comece com um limite baixo: a interface mostra, ao fim da coleta, a cobertura obtida por campo e a lista de atributos e valores que o vocabulário não reconheceu. Acrescente as entradas pertinentes pelos editores da própria página e execute de novo.

## Vocabulários

| Arquivo | Conteúdo |
| --- | --- |
| `vocabulario.csv` | `campo_final,sinonimo` — une nomes de atributo divergentes em um campo |
| `valores.csv` | `campo,valor_final,sinonimo` — padroniza os valores dentro de um campo |
| `indicios.csv` | `campo,padrao` — expressões regulares procuradas nos campos textuais |

São lidos a cada execução. Editá-los não exige mexer no código.

## Limitações

As métricas de desempenho publicadas no TCC foram obtidas contra um conjunto de 187 amostras de *Astyanax mexicanus* curadas manualmente, provenientes de 11 BioProjects. Cada campo é testável em uma a três submissões independentes, de modo que os valores descrevem o comportamento do método diante de um conjunto restrito de convenções, não sua generalidade.

O vocabulário em sua versão final foi construído com acesso a esse mesmo conjunto, o que caracteriza desempenho em ajuste e não em generalização.

Os rótulos dos campos finais são locais, em português, e não estão mapeados a identificadores de ontologia.

## Licença

MIT. Ver `LICENSE`.

## Como citar

Ver `CITATION.cff`.
