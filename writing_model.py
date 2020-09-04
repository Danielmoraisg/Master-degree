#escrever os metabolitepools
import os
directory = os.getcwd()
folder = list(next(os.walk(directory)))[1]
import pandas as pd
relacao = {}
df = pd.DataFrame()
lista = []
for pasta in folder:
    new_folders = list(next(os.walk("%s\\%s"%(directory,pasta))))
    for file in new_folders[2]:
        if ".txt" in file and "ncomponentes_" in file:
            data = pd.read_csv("%s\\%s"%(new_folders[0],file))
            relacao.update(dict(zip(data["id"],data["carbon_number"])))
            df = df.append(data,ignore_index = True)
            df = df.drop_duplicates()
#pegar a versao e somar 1
try:
    arquivo = open("model.txt","r")
    for linha in arquivo.read().split("\n"):
        if "<version>" in linha:
            version = float(linha[11:14])+1
            break
except:
    version = 1
arquivo = open("model.txt","w")
arquivo.write('<?xml version="1.0" encoding="utf-8"?>\n<fluxml xmlns="http://www.13cflux.net/fluxml"> \n\t<info> \n\t\t<version>%.1f</version> \n\t\t<comment/> \n\t</info> \n\t<reactionnetwork> \n\t\t<metabolitepools>\n'%version)
cofactors = ["H","NADP","NADPH","NAD","NADH","ATP","ADP","QUINOL","QUINONA","UDP","UTP","COA"]
relacao['ARG'] = 6
for compound in relacao:
    if relacao[compound] != 0 and compound not in cofactors:
        arquivo.write('\t\t\t<pool id = "%s" atoms = "%s"/>\n'%(compound,relacao[compound]))
arquivo.write("\t\t</metabolitepools>\n")
#puxar dataframe com os dados de reversibilidade
df_rev = pd.read_csv("C:/Users/danie/Desktop/mestrado/Model/df_rev.csv")
dict_rev = {0.0:"true", 1.0:"false"}
dict_rev_compounds = dict(zip(df_rev["nome_generico"],df_rev["reversibilidade"]))#escrever as reactions
#escrever as reactions
for pasta in folder:
    new_folders = list(next(os.walk("%s\\%s"%(directory,pasta))))
    for file in new_folders[2]:
        if "nReacoes_" in file:
            data = pd.read_csv("%s\\%s"%(new_folders[0],file))
            for i in range(len(data)):
                dados = list(data.iloc[i])
                lista.append(dados[0])
                try:
                    arquivo.write('\t\t<reaction id="%s" bidirectional = "%s">\n'%(dados[0], dict_rev[dict_rev_compounds[dados[0]]]))
                except KeyError:
                    arquivo.write('\t\t<reaction id="%s" bidirectional = "false">\n'%(dados[0]))
                if ";" in str(dados[2]):
                    reduct_id = str(dados[2]).split(";")
                    reduct_config = str(dados[1]).split(";")
                    for count in range(len(reduct_id)):
                        if reduct_id[count] not in cofactors:
                            arquivo.write('\t\t\t<reduct cfg="%s" id="%s"/>\n'%(reduct_config[count],reduct_id[count]))
                else:
                    reduct_id = dados[2]
                    reduct_config = dados[1]
                    arquivo.write('\t\t\t<reduct cfg="%s" id="%s"/>\n'%(reduct_config,reduct_id))
                if ";" in str(dados[4]):
                    rproduct_id = str(dados[4]).split(";")
                    rproduct_config = str(dados[3]).split(";")
                    for count in range(len(rproduct_id)):
                        if rproduct_id[count] not in cofactors:
                            arquivo.write('\t\t\t<rproduct cfg="%s" id="%s"/>\n'%(rproduct_config[count],rproduct_id[count]))
                else:
                    rproduct_id = dados[4]
                    rproduct_config = dados[3]
                    arquivo.write('\t\t\t<rproduct cfg="%s" id="%s"/>\n'%(rproduct_config,rproduct_id))
                arquivo.write('\t\t</reaction>\n')
#escrever os efluxos udp, precursores de ac hia, co2, trehalose
efflux = ["CO2","UDPNACETYLG","UDPAC_GLU","ET","GLY","ACET","VAL","ALA","ASP","ASN",'THRE','MET','GLUT','GLI','ARG','SER','GLN','PHE','TYR','LEU','ILE','LYS','HIS','CYS','PRO','TRP','ACOA','MetTHF']
alphabet = []
for letter in range(97,123):
    alphabet.append(chr(letter))
count = 0
for compound in efflux:
    arquivo.write('\t\t<reaction id="Efflux%s">\n'%(count))
    count+=1
    arquivo.write('\t\t\t<reduct cfg="%s" id="%s"/>\n'%("".join(alphabet[0:int(df.loc[df["id"] == compound]["carbon_number"])]),compound))
    arquivo.write('\t\t</reaction>\n')
arquivo.write('\t</reactionnetwork>\n')
#escrevendo o input
arquivo.write('\t<configuration name="default" stationary="true">\n\t\t<input pool="GLC1" type="isotopomer">\n\t\t\t<label cfg="110111">0.0016</label>\n\t\t\t<label cfg="011111">0.0016</label>\n\t\t\t<label cfg="101111">0.0016</label>\n\t\t\t<label cfg="111101">0.0016</label>\n\t\t\t<label cfg="111111">0.9904</label>\n\t\t\t<label cfg="111110">0.0016</label>\n\t\t\t<label cfg="111011">0.0016</label>\n\t\t</input>\n\t\t<input pool="GLC0" type="isotopomer">\n\t\t\t<label cfg="000100">0.011</label>\n\t\t\t<label cfg="000000">0.934</label>\n\t\t\t<label cfg="000001">0.011</label>\n\t\t\t<label cfg="000010">0.011</label>\n\t\t\t<label cfg="001000">0.011</label>\n\t\t\t<label cfg="010000">0.011</label>\n\t\t\t<label cfg="100000">0.011</label>\n\t\t</input>')
#corrigindo simbolos do modelo
arquivo = open("model.txt", "r")
data = arquivo.read()
arquivo.close()
proibidos = ["-","_","Glycolysis0G"]
for i in proibidos:
    if i == "Glycolysis0G":
        data = data.replace(i,"2PG")
    else:
        data = data.replace(i,"")
arquivo = open("model.txt", "w")
arquivo.write(data)

#escrevendo constraints
#importando os documentos
arquivo.write('\n\t<constraints>')
arquivo.write('\n\t\t<net>')
arquivo.write('\n\t\t\t<textual>')
doc_constraints = open("C:/Users/danie/Desktop/mestrado/Model/constraints/Constraints.txt", "r")
df_taxas = pd.read_csv("C:/Users/danie/Desktop/mestrado/Model/constraints/taxas.txt", sep = ",")
#criando onde as taxas entram em relação às reações
dict_rendimento = dict(zip(df_taxas.compounds,df_taxas.produtividade_especifica))
rel_compound_reaction = {"Glicose":"Glycolysis4", "Acetato":"acetate0", "Glicerol":"Glycerol0","Etanol":"etanol1", "biomassa":"BM0"}
#colocando constraints no documento utilizando as constraints calculadas(taxas) e escritas (proporções)
constraints = []
for compound in dict_rendimento:
    constraints.append("%s=%s"%(rel_compound_reaction[compound],dict_rendimento[compound]))
for constraint in doc_constraints.read().split("\n")[1:]:
    constraints.append(constraint)
arquivo.write(";".join(constraints))
arquivo.write('</textual>')
arquivo.write('\n\t\t</net>')
arquivo.write('\n\t</constraints>')
doc_constraints.close()
#importando file com o output do isocor
import pandas as pd
os.chdir(r'C:\Users\danie\Desktop\mestrado\CG-MS\Isocor')
df = pd.read_csv("data_full_res.tsv", sep = "\t")
df = df.append(pd.read_csv("data_full_cromato_res.tsv", sep = "\t"),ignore_index = True)
#subset do dataframe com as colubas metabolite isotopologue e isotopologue fraction sem os NANs
isocor_df = df[['metabolite','isotopologue', 'isotopologue_fraction']].dropna()
#checar as posicoes de marcacao isotopica na tabela tbl_fragmantation
os.chdir(r"C:\Users\danie\Desktop\mestrado\CG-MS")
df = pd.read_csv("Tbl_fragmentation.csv",sep = ";", skiprows = 3)
df_fragmentation = subdf = df[['Nameg\n\n(MEOX, methoxyamination)\n(TMS, trimethylsilylation)\n','Mass to Charge Ratio of \nMonoisotopic \n12C Mass Isotopomers','C-Positions within \nMetabolite']]
#Trocar o nome das colunas
df_fragmentation.columns = ['composto','mz','C_positions']
#relacionar fragmento com conjundo de carbonos labelled
frag = []
cont = 0
for i in df_fragmentation["composto"]:
    fragment = '-'.join(list(map(str,df_fragmentation.loc[cont,['composto','mz']])))
    frag.append(fragment)
    cont+=1
rel = []
for i in frag:
    for j in isocor_df["metabolite"].drop_duplicates():
        if (j == i[0:j.rfind("-")]+i[i.rfind("-"):]):
            rel.append(i)
#tratar glutamina:
for i in rel:
    if 'Glutamine' in i :
        rel.remove(i)
#colocando measurements
os.chdir(r'C:\Users\danie\Desktop\mestrado\Model')
arquivo.write("\n\t\t<measurement>")
arquivo.write("\n\t\t\t<model>")
arquivo.write('\n\t\t\t\t<labelingmeasurement>')
group = 0
index = 0
rel_sigla = {'Serine (2TMS)':'Ser','Threonine (3TMS)':'Thre',"Citric acid (4TMS)e":"Cit",'Malic acid (3TMS)':'Mal', 'Alanine (2TMS)':'Ala','Glutamic acid (3TMS)': 'Glut', "Trehalose, alpha,alpha'- (8TMS)": 'Tre', 'Tyrosine (3TMS)': 'Tyr', 'Valine (2TMS)': 'Val', 'Succinic acid (2TMS)': 'Suc','Fumaric acid (2TMS)':"Fuma", 'Glucose (1MEOX) (5TMS)': 'Glu','Leucine (2TMS)':'Leu','Glycine (2TMS)':'Gli'}
rel_frag_c = dict(zip(frag,df_fragmentation["C_positions"]))
ms_groups = {}
#escrevendo o campo
import math
for comp in rel:
    try:
        compound = rel_sigla[comp[0:comp.rfind('-')]]+comp[comp.rfind('-'):]
    except KeyError:
        continue
    if  compound[0:compound.find('-')].upper() in relacao.keys():
        group+=1
        try :
            math.isnan(rel_frag_c[comp])
        except TypeError:
            ms_groups['"ms_group_%s"'%group]=compound
            arquivo.write('\n\t\t\t\t\t<group id="ms_group_%s" scale="auto">\n\t\t\t\t\t\t<textual>%s[%s]#M%s</textual>\n\t\t\t\t\t</group>'%(group,compound[0:compound.find('-')].upper(),rel_frag_c[comp],",".join(map(str,list(range(isocor_df.loc[isocor_df.where(isocor_df == compound).last_valid_index(),'isotopologue']+1))))))
    else:
        print('O composto"%s" nao esta no modelo!'%compound)
        group+=1
        ms_groups['"ms_group_%s"'%group]=compound
        arquivo.write('\n\t\t\t\t\t<group id="ms_group_%s" scale="auto">\n\t\t\t\t\t\t<textual>%s[%s]#M%s</textual>\n\t\t\t\t\t</group>'%(group,compound[0:compound.find('-')].upper(),rel_frag_c[comp],",".join(map(str,list(range(isocor_df.loc[isocor_df.where(isocor_df == compound).last_valid_index(),'isotopologue']+1))))))

    index+=1
arquivo.write("\n\t\t\t\t</labelingmeasurement>")
arquivo.write("\n\t\t\t</model>")
#escrever dados do isocor
arquivo.write("\n\t\t\t<data>")
for group in ms_groups:
    weight = 0
    subdf = isocor_df[isocor_df.metabolite == ms_groups[group]]
    for compound in subdf["isotopologue_fraction"]:
        try:
            arquivo.write('\n\t\t\t\t<datum id=%s stddev = "0.020" weight="%i">%.3f</datum>'%(group,weight,float(subdf[subdf.isotopologue == weight]['isotopologue_fraction'].tail(1))))
        except TypeError:
            break
        weight+=1
arquivo.write("\n\t\t\t</data>")
arquivo.write("\n\t\t</measurement>")
arquivo.write('\n\t</configuration>')
#terminando
arquivo.write("\n</fluxml>")
arquivo.close()
#arrumar os nomes e colocar sebrepor AA no metabolismo central de carbono:
arquivo = open("model.txt", "r")
data = arquivo.read()
arquivo.close()
data = data.replace('3PG','TPG')
data = data.replace('2OXOMIT','DOXOMIT')
arquivo = open("model.txt", "w")
arquivo.write(data)
arquivo.close()

