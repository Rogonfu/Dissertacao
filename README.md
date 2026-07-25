# Dissertacao

Para compreender adequadamente o conteúdo disponibilizado neste repositório e seu contexto, recomenda-se a leitura prévia da dissertação.

O arquivo "pre_processamento" realiza o pré-processamento dos conjuntos de dados obtidos no repositório GEO (https://www.ncbi.nlm.nih.gov/geo/). Para utilizá-lo, é necessário realizar o download dos conjuntos de dados por meio de seus respectivos identificadores (IDs). Informações detalhadas sobre essa etapa encontram-se no capítulo "Desenvolvimento do Trabalho" da dissertação. Após a execução desse script, os conjuntos de dados pré-processados foram armazenados para evitar a necessidade de repetir essa etapa em futuras execuções. Os dez conjuntos de dados pré-processados utilizados nos experimentos estão disponíveis na pasta "DATA".

O arquivo "GENES.zip" contém todas as listas de genes extraídas da base de dados DisGeNET, correspondentes às patologias analisadas nos conjuntos de dados da GEO, bem como as informações dos genes e de suas respectivas famílias obtidas na base HGNC.

O arquivo "Algoritmo.ipynb" contém o carregamento dos dados presentes no arquivo "GENES.zip", bem como dos conjuntos de dados pré-processados localizados na pasta "DATA", seguido da execução das etapas do framework GSM e implementação dos algoritmos de seleção de atributos baseado em filtro.

O arquivo "Algoritmo-GSMC.ipynb" contém todas as etapas presentes no arquivo "Algoritmo.ipynb" com a exeção dos algoritmos baseado em filtro, acrescendo apenas a nova etapa de combinação de subconjuntos incorporada ao framework GSM.



