# Dissertacao

Para compreender todo o material aqui presente e o seu contexto, deve-se ler a dissertação primeiramente.

O arquivo "pre_processamento" realiza o pré-pocessamento dos conjuntos de dados obtidos do repositório GEO (https://www.ncbi.nlm.nih.gov/geo/). Dentro do repositório deve-se buscar os conjuntos de dados pelos seus respectivos IDs. Para mais informações, deve-se ler o capítulo "Desenvolvimento do trabalho" na dissertação. Após a execução deste script, foram guardados os conjuntos de dados pré-processados para não precisar executar o código novamente todas as vezes. Todos os 10 conjuntos de dados pré-processados estão dentro da pasta "DATA" e foram utilizados nos algoritmos desenvolvidos.

O arquivo "GENES.zip" contém todas as listas de genes extraidas da base de dados DisGeNet para cada patologia dos conjuntos de dados da GEO e também todas as informações dos genes e suas famílias extraidas da base HGNC.

O arquivo "Algoritmo.ipynb" contém o carregamento dos conjuntos zipados no arquivo "GENES.zip", junto do carregamento dos conjuntos de dados pré-processados da pasta DATA e as subsequentes etapas do framework GSM.

O arquivo "Algoritmo-GSMC" contém as mesmas etapas do arquivo "Algoritmo.ipynb", modificando apenas a nova etapa de combinação adicionada ao framework.



