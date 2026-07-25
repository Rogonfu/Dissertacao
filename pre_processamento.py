import pandas as pd
import numpy as np
from imblearn.under_sampling import RandomUnderSampler
from collections import defaultdict
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import ShuffleSplit
import math
from sklearn.model_selection import cross_val_score

d1 = pd.read_csv("GENES\\family.csv", sep=",", on_bad_lines='warn')
d2 = pd.read_csv("GENES\\gene_has_family.csv", sep=",", on_bad_lines='warn')
d3 = pd.read_csv("GENES\\genes.tsv", sep="\t", on_bad_lines='warn')

df1 = pd.DataFrame()

def Conjunto():
    caminho = input("Digite o número do conjunto de dados a ser carregado: \n\n1- Ulcerative Colitis\n2- Glioma\n3- Metastatic prostate cancer (HG-U95C)\n4- Metastatic prostate cancer (HG-U95A)\n5- Lung cancer\n6- Lung adenocarcinoma\n7- Leukemia\n8- Pulmonary hypertension\n9- Non-small cell lung carcinoma\n10- Colorectal cancer \n\n")

    global df1

    if caminho == '1':       # ulcerative colitis
        df = pd.read_csv('DATA\\GDS3268.soft', sep="\t",header=None, on_bad_lines='warn')

        #Genes
        df1 = pd.read_csv("GENES\\colitis.tsv", sep="\t", on_bad_lines='warn')

        # Atribuição de números as classes
        ot = ['GSM282855','GSM282856','GSM282857','GSM282858','GSM282859','GSM282860','GSM282861','GSM282862','GSM282863','GSM282864','GSM282865','GSM282866','GSM282867','GSM282868','GSM282869','GSM282870','GSM282871','GSM282872','GSM282873','GSM282874','GSM282875','GSM282876','GSM282877','GSM282878','GSM282879','GSM282880','GSM282881','GSM282882','GSM282883','GSM282884','GSM282885','GSM282886','GSM282887','GSM282888','GSM282889','GSM282890','GSM282891','GSM282892','GSM282893','GSM282894','GSM282895','GSM282896','GSM282897','GSM282898','GSM282899','GSM282900','GSM282901','GSM282902','GSM282903','GSM282904','GSM282905','GSM282906','GSM282907','GSM282908','GSM282909','GSM282910','GSM282911','GSM282912','GSM282913','GSM282914','GSM282915','GSM282916','GSM282917','GSM282918','GSM282919','GSM282920','GSM282921','GSM282922','GSM282923','GSM282924','GSM282925','GSM282926','GSM282927']

    elif caminho == '2':         # Glioma
        df = pd.read_csv("DATA\\GDS1962.soft", sep="\t",header=None, on_bad_lines='warn')

        #Genes
        df1 = pd.read_csv("GENES\\Glioma.tsv", sep="\t", on_bad_lines='warn')

       # Atribuição de números as classes
        ot = ['GSM97800','GSM97803','GSM97804','GSM97805','GSM97807','GSM97809','GSM97811','GSM97812','GSM97816','GSM97817','GSM97820','GSM97825','GSM97827','GSM97828','GSM97833','GSM97834','GSM97840','GSM97846','GSM97848','GSM97849','GSM97850','GSM97853','GSM97855']

    elif caminho == '3':      # Metastatic prostate cancer (HG-U95C)
        df = pd.read_csv("DATA\\GDS2547.soft", sep="\t", header=None, on_bad_lines='warn')

        # Genes
        df1 = pd.read_csv("GENES\\Prostate.tsv", sep="\t", on_bad_lines='warn')

        # Atribuição de números as classes
        ot = ['GSM152839','GSM152840','GSM152841','GSM152842','GSM152843','GSM152844','GSM152845','GSM152846','GSM152847','GSM152848','GSM152849','GSM152850','GSM152851','GSM152852','GSM152853','GSM152854','GSM152855','GSM153238','GSM153239','GSM153240','GSM153241','GSM153242','GSM153243','GSM153244','GSM153245','GSM153246','GSM153247','GSM153248','GSM153249','GSM153250','GSM153251','GSM153252','GSM153253','GSM153254','GSM153255','GSM153256','GSM153257','GSM153258','GSM153259','GSM153260','GSM153261','GSM153262','GSM153263','GSM153264','GSM153265','GSM153266','GSM153267','GSM153268','GSM153269','GSM153270','GSM153271','GSM153272','GSM153273','GSM153274','GSM153275','GSM153276','GSM153277','GSM153278','GSM153279','GSM153280','GSM153281','GSM153282','GSM153283','GSM153284','GSM153285','GSM153286','GSM153287','GSM153288','GSM153289','GSM153290','GSM153291','GSM153292','GSM153293','GSM153294','GSM153295']

    elif caminho == '4':       # Metastatic prostate cancer (HG-U95A)
        df = pd.read_csv("DATA\\GDS2545.soft", sep="\t",header=None, on_bad_lines='warn')

        # Genes
        df1 = pd.read_csv("GENES\\Prostate.tsv", sep="\t", on_bad_lines='warn')

        # Atribuição de números as classes
        ot = ['GSM152804','GSM152805','GSM152806','GSM152807','GSM152808','GSM152809','GSM152810','GSM152811','GSM152812','GSM152813','GSM152814','GSM152815','GSM152816','GSM152817','GSM152818','GSM152819','GSM152820','GSM152821','GSM153115','GSM153116','GSM153117','GSM153118','GSM153119','GSM153120','GSM153121','GSM153122','GSM153123','GSM153124','GSM153125','GSM153126','GSM153127','GSM153128','GSM153129','GSM153130','GSM153131','GSM153132','GSM153133','GSM153134','GSM153135','GSM153136','GSM153137','GSM153138','GSM153139','GSM153140','GSM153141','GSM153142','GSM153143','GSM153144','GSM153145','GSM153146','GSM153147','GSM153148','GSM153149','GSM153150','GSM153151','GSM153152','GSM153153','GSM153154','GSM153155','GSM153156','GSM153157','GSM153158','GSM153159','GSM153160','GSM153161','GSM153162','GSM153163','GSM153164','GSM153165','GSM153166','GSM153167','GSM153168','GSM153169','GSM153170','GSM153171','GSM153172','GSM153173','GSM153174','GSM153175','GSM153176','GSM153177']

    elif caminho == '5':   # lung cancer
        df = pd.read_csv("DATA\\GDS2771.soft", sep="\t",header=None, on_bad_lines='warn')

        # Genes
        df1 = pd.read_csv("GENES\\LungCarcinoma.tsv", sep="\t", on_bad_lines='warn')

        ot = ['GSM93997','GSM94077','GSM94078','GSM94079','GSM94080','GSM94081','GSM94082','GSM94083','GSM94084','GSM94085','GSM94086','GSM94087','GSM94088','GSM94089','GSM94090','GSM94091','GSM94092','GSM94093','GSM94094','GSM94095','GSM94096','GSM94097','GSM94098','GSM94099','GSM94100','GSM94101','GSM94102','GSM94103','GSM94104','GSM94105','GSM94106','GSM94107','GSM94108','GSM94109','GSM94110','GSM94111','GSM94112','GSM94113','GSM94114','GSM94115','GSM94116','GSM94117','GSM94118','GSM94119','GSM94120','GSM94121','GSM94122','GSM94123','GSM94124','GSM94125','GSM94126','GSM94127','GSM94128','GSM94129','GSM94130','GSM94131','GSM94132','GSM94133','GSM94134','GSM94135','GSM94136','GSM94137','GSM94138','GSM94139','GSM94140','GSM94141','GSM94142','GSM94143','GSM94144','GSM94145','GSM94146','GSM94147','GSM94148','GSM98785','GSM98786','GSM98787','GSM98788','GSM98789','GSM98790','GSM98791','GSM98792','GSM98793','GSM98794','GSM98795','GSM98796','GSM98797','GSM98798','GSM98799','GSM98800','GSM98801']

    elif caminho == '6':       # Lung adenocarcinoma
        df = pd.read_csv("DATA\\GDS3257.soft", sep="\t",header=None, on_bad_lines='warn')

        # Genes
        df1 = pd.read_csv("GENES\\AdenocarcinomaLung.tsv", sep="\t", on_bad_lines='warn')

        ot = ['GSM254626','GSM254628','GSM254632','GSM254634','GSM254635','GSM254638','GSM254640','GSM254643','GSM254644','GSM254646','GSM254649','GSM254651','GSM254653','GSM254655','GSM254658','GSM254660','GSM254662','GSM254665','GSM254667','GSM254669','GSM254671','GSM254673','GSM254676','GSM254677','GSM254679','GSM254681','GSM254683','GSM254685','GSM254689','GSM254691','GSM254693','GSM254695','GSM254699','GSM254702','GSM254703','GSM254706','GSM254708','GSM254710','GSM254711','GSM254712','GSM254713','GSM254715','GSM254717','GSM254719','GSM254723','GSM254725','GSM254727','GSM254730','GSM254731']

    elif caminho == '7':       # Leukemia
        df = pd.read_csv("DATA\\GDS4206.soft", sep="\t",header=None, on_bad_lines='warn')

        # Genes
        df1 = pd.read_csv("GENES\\leukemia.tsv", sep="\t", on_bad_lines='warn')

        ot = ['GSM342328','GSM342418','GSM342419','GSM342420','GSM342422','GSM342423','GSM342424','GSM342425','GSM342427','GSM342428','GSM342429','GSM342430','GSM342431','GSM342432','GSM342433','GSM342434','GSM342435','GSM342436','GSM342330','GSM342438','GSM342439','GSM342440','GSM342442','GSM342443','GSM342445','GSM342446','GSM342447','GSM342448','GSM342449','GSM342450','GSM342452','GSM342453','GSM342455','GSM342456','GSM342457','GSM342332','GSM342458','GSM342459','GSM342460','GSM342461','GSM342463','GSM342464','GSM342465','GSM342466','GSM342333','GSM342468','GSM342469','GSM342470','GSM342472','GSM342473','GSM342474','GSM342475','GSM342476','GSM342478','GSM342479','GSM342480','GSM342481','GSM342482','GSM342483','GSM342485','GSM342486','GSM342487','GSM342488','GSM342489','GSM342490','GSM342492','GSM342493','GSM342494','GSM342497','GSM342336','GSM342498','GSM342499','GSM342500','GSM342501','GSM342502','GSM342503','GSM342504','GSM342505','GSM342506','GSM342507','GSM342337','GSM342508','GSM342509','GSM342510','GSM342511','GSM342512','GSM342513','GSM342514','GSM342515','GSM342320','GSM342339','GSM342340','GSM342341','GSM342342','GSM342343','GSM342344','GSM342345','GSM342346','GSM342321','GSM342348','GSM342349','GSM342350','GSM342351','GSM342352','GSM342353','GSM342355','GSM342358','GSM342359','GSM342360','GSM342361','GSM342362','GSM342363','GSM342364','GSM342365','GSM342367','GSM342368','GSM342369','GSM342370','GSM342372','GSM342373','GSM342374','GSM342375','GSM342376','GSM342377','GSM342324','GSM342378','GSM342379','GSM342381','GSM342383','GSM342384','GSM342386','GSM342387','GSM342388','GSM342389','GSM342390','GSM342391','GSM342392','GSM342393','GSM342394','GSM342395','GSM342396','GSM342397','GSM342398','GSM342400','GSM342401','GSM342403','GSM342405','GSM342406','GSM342407','GSM342410','GSM342411','GSM342412','GSM342413','GSM342414','GSM342415','GSM342416','GSM342417']

    elif caminho == '8':      # Pulmonary hypertension
        df = pd.read_csv("DATA\\GDS5499.soft", sep="\t",header=None, on_bad_lines='warn')

        # Genes
        df1 = pd.read_csv("GENES\\pulmonary.tsv", sep="\t", on_bad_lines='warn')

        ot = ['GSM827665','GSM827666','GSM827667','GSM827668','GSM827669','GSM827670','GSM827671','GSM827672','GSM827673','GSM827674','GSM827675','GSM827676','GSM827677','GSM827678','GSM827679','GSM827680','GSM827681','GSM827682','GSM827683','GSM827684','GSM827685','GSM827686','GSM827687','GSM827688','GSM827689','GSM827690','GSM827691','GSM827692','GSM827693','GSM827694','GSM827695','GSM827696','GSM827697','GSM827698','GSM827699','GSM827700','GSM827701','GSM827702','GSM827703','GSM827704','GSM827705']

    elif caminho == '9':       # Non-small cell lung carcinoma
        df = pd.read_csv("DATA\\GDS3837.soft", sep="\t",header=None, on_bad_lines='warn')

        # Genes
        df1 = pd.read_csv("GENES\\non-small.tsv", sep="\t", on_bad_lines='warn')

        ot = ['GSM494616','GSM494617','GSM494618','GSM494619','GSM494620','GSM494621','GSM494622','GSM494623','GSM494624','GSM494625','GSM494626','GSM494627','GSM494628','GSM494629','GSM494630','GSM494631','GSM494632','GSM494633','GSM494634','GSM494635','GSM494636','GSM494637','GSM494638','GSM494639','GSM494640','GSM494641','GSM494642','GSM494643','GSM494644','GSM494645','GSM494646','GSM494647','GSM494648','GSM494649','GSM494650','GSM494651','GSM494652','GSM494653','GSM494654','GSM494655','GSM494656','GSM494657','GSM494658','GSM494659','GSM494660','GSM494661','GSM494662','GSM494663','GSM494664','GSM494665','GSM494666','GSM494667','GSM494668','GSM494669','GSM494670','GSM494671','GSM494672','GSM494673','GSM494674','GSM494675']

    elif caminho == '10':     # Colorectal cancer
        df0 = pd.read_csv("DATA\\GDS4516.soft", sep="\t",header=None, on_bad_lines='warn')
        df = pd.read_csv("DATA\\GDS4718.soft", sep="\t",header=None, on_bad_lines='warn')
        #df = pd.concat([df0, df], axis=1, join='inner')
        df = pd.merge(df0, df, on=[0, 1], how="inner")

        # Genes
        df1 = pd.read_csv("GENES\\colorectal.tsv", sep="\t", on_bad_lines='warn')

        ot = ['GSM549121','GSM549102','GSM549104','GSM549108','GSM549119','GSM549133','GSM549139','GSM549099','GSM549109','GSM549110','GSM549114','GSM549122','GSM549134','GSM549136','GSM549140','GSM549111','GSM549113','GSM549132','GSM549137','GSM549142','GSM549100','GSM549107','GSM549115','GSM549116','GSM549120','GSM549131','GSM549118','GSM549129','GSM549123','GSM549124','GSM549126','GSM549128','GSM549103','GSM549117','GSM549138','GSM549141','GSM549130','GSM549101','GSM549105','GSM549106','GSM549112','GSM549125','GSM549127','GSM549135']

    else:
        return 0


    print('Conjunto de dados original:\n')
    print(df)

    #pré-processamento
    new = df.dropna(subset=[1]).drop(0, axis=1).set_index(1).T.dropna(axis=1)

    new['IDENTIFIER'] = new['IDENTIFIER'].astype(str).str.strip()

    # atribui 0 quando o IDENTIFIER está na lista ot, caso contrário 1
    new['IDENTIFIER'] = np.where(new['IDENTIFIER'].isin(ot), 0, 1)

    clas = new['IDENTIFIER']
    conjunto = new.drop('IDENTIFIER', axis=1).astype('float')
    clas = clas.astype('int')
    conjunto['Classe'] = clas

    #under sampling
    rus = RandomUnderSampler(random_state=42)
    X_res, y_res = rus.fit_resample(conjunto.drop('Classe', axis=1), clas)
    conjunto = X_res
    conjunto['Classe'] = y_res
    clas = conjunto['Classe']

    # Atualizar genes do conjunto
    d4 = d3.loc[:, ['hgnc_id', 'alias_symbol']]
    d4 = d4.dropna()

    dic = defaultdict(list)

    for gene in conjunto.columns:
        for row in d4.itertuples():
            ids = row.hgnc_id

            subgenes = row.alias_symbol
            subgenes = subgenes.split('|')

            for i in subgenes:
                if gene == i:
                    dic[ids].append(gene)

    for chave, valor in dic.items():
        c = set(valor)
        var = d3.loc[d3['hgnc_id'] == chave]
        if len(c) == 1:
            conjunto.rename(columns={valor[0]: var.iat[0,1]}, inplace=True)
        else:
            for i in c:
                conjunto.rename(columns={i: var.iat[0,1]}, inplace=True)

    conjunto = conjunto.sample(frac = 1, random_state=42)
    conjunto = conjunto.reset_index(drop=True)

    print("Conjunto de dados pré-processado:\n")
    print(conjunto)

    return conjunto