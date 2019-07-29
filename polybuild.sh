#!/bin/bash

#=========================================================================================================
# Programa que escreve um arquivo de input, a ser lido pelo Playmol, com as coordenadas
# de um sistema polimérico polidisperso.

# Desenvolvido por Tiago S. M. Lemos - PEQ/COPPE/UFRJ
# Rio de Janeiro - RJ, 1º de outubro de 2016
#
# Edições - incluída a possibilidade de novos componentes (misturas)
#   Versão 0.1 - 05/10/2016
#   o arquivo com a distribuição tem a seguinte estrutura:
#   nome do componente tamanho de cadeia (1 para moléculas não oligoméricas quantidade

#   Versão 0.2 - 06/10/2016
#   Escrita do arquivo Lammps, para a simulação do sistema gerado
#   Possibilidade de escrita de copolímeros dibloco

#   Versão 1.0 - 16/10/2016
#   Substituição do tipo de input: ao invés de somente a microestrura, o input agora 
#   contém as informações sobre todos os parâmetros relevantes à simulação

#========================================================================================================


# Nome dos arquivos de entrada (Distribuição de Tamanhos de Cadeias - CLD) e de saída (input Playmol):

for ParameterFile in $*
do

echo "arquivo de configuração ${ParameterFile} sendo processado pelo polybuild"


# Crie um diretório para incluir os resultados e mova para lá os arquivos necessários a execução do programa:

prefixo=$(echo "${ParameterFile}" | cut -d. -f 1)
if ! test -d "${prefixo}_results"; then
  mkdir "${prefixo}_results"
fi
cp ${ParameterFile} "${prefixo}_results"
cd "${prefixo}_results"
cp ${ParameterFile} parameterfile.backup


# Remova os espaços em branco repetidos e as linhas em branco do arquivo de estrutura:

sed 's/[ ]\+/ /g ; /^ *\?$/d' parameterfile.backup > ${ParameterFile}




# Defina os nomes dos arquivos de saída:

InputPlaymolFile="${prefixo}.mol"                  # nome do arquivo de input a ser lido pelo Playmol
xyz_file="${prefixo}.xyz"                          # nome do arquivo xyz a ser gerado
vmd_read_file="${prefixo}.lammpstrj"               # nome do arquivo passível de visualização por VMD
lammps_data_file="data.${prefixo}.lmp"             # nome do arquivo de dados gerado para a simulação Lammps 


# Definição das variáveis do sistema a ser simulado:

k_spring=$(            grep -E "^k_spring "               ${ParameterFile}  | cut -d " " -f 3)
xeq_spring=$(          grep -E "^xeq_spring "             ${ParameterFile}  | cut -d " " -f 3)


# Definição dos parâmetros Playmol:                                                                  
ro=$(                  grep -E "^ro "                   ${ParameterFile}  | cut -d " " -f 3)
r_bond=$(              grep -E "^r_bond "               ${ParameterFile}  | cut -d " " -f 3)
seed=$(                grep -E "^seed "                 ${ParameterFile}  | cut -d " " -f 3)
tol=$(                 grep -E "^tol "                  ${ParameterFile}  | cut -d " " -f 3)
retry_value=$(         grep -E "^retry_value "          ${ParameterFile}  | cut -d " " -f 3)
                                                                                                               

# Reconhecimento dos tipos de átomos e criação de uma lista com seus nomes:
init_micro=$(($(    grep -E -n "Microestrutura"           ${ParameterFile}  | cut -d ":" -f 1)+1))
final_micro=$(($(   grep -E -n "Lista de componentes"     ${ParameterFile}  | cut -d ":" -f 1)-1)) 
init_compmass=$(($( grep -E -n "Lista de componentes"     ${ParameterFile}  | cut -d ":" -f 1)+1)) 
final_compmass=$(($(grep -E -n "Coeficiente de interação" ${ParameterFile}  | cut -d ":" -f 1)-1)) 

sed -n ${init_micro},${final_micro}p        ${ParameterFile} > structure.tmp
sed -n ${init_compmass},${final_compmass}p  ${ParameterFile} > components.tmp

# Escrita de informações iniciais dos inputs Playmol:

echo -e "**********Escrevendo input para o Playmol********** \n. \n.. \n... "
echo "# Arquivo gerado a partir do script 'inp_play_homopdi'
# data de criação: $(date) 

  "  > ${InputPlaymolFile}

while read id tipo massa; do
  echo "atom_type       ${tipo} 
mass            ${tipo}       ${massa}"  >> ${InputPlaymolFile}
done < components.tmp
rm components.tmp

echo "define          ro  as ${ro}            # numeric density
bond_type       *   *   ${k_spring} ${xeq_spring}
  
"  >> ${InputPlaymolFile}


# Leitura do arquivo com as informações da CLD e escrita do arquivo xyz:

atom_count=1
molecule_count=1
> file1.tmp
> ${xyz_file}
while read line; do

  echo "# molécula  ${molecule_count} " >> ${InputPlaymolFile}
  n_words="$(echo "${line}" | wc -w)"

  quantidade="$(echo "${line}" | cut -d " " -f ${n_words})"
  echo "# N${molecule_count} =  ${quantidade}" >> ${InputPlaymolFile}
  echo "packmol    copy ${molecule_count}  ${quantidade} " >> file1.tmp
  molecule_count=$((${molecule_count}+1))

  atom_init=${atom_count}
  x_atom=0.000 

  for n in $(seq 1 2 $((${n_words}-1))) ; do
    tipo="$(echo "${line}" | cut -d " " -f ${n})"
    tamanho_bloco="$(echo "${line}" | cut -d " " -f $((${n}+1)))"

    for n2 in $(seq 1 ${tamanho_bloco}); do 
      echo "atom     ${atom_count}    ${tipo}" >> ${InputPlaymolFile}
      echo "${atom_count}      ${x_atom}    0.000    0.000" >> ${xyz_file}
      # variavel=$(echo "expressão" | bc)            # expressão para uso do bc
      x_atom=$(echo " scale=3; ${x_atom}+${r_bond}" | bc)
      atom_count=$((${atom_count}+1))
    done
  done
  atom_final=$((${atom_count}-1))

  echo "" >> ${InputPlaymolFile}
  for n in $(seq ${atom_init} $(( ${atom_final}-1)) ); do
    if  ((${atom_init}!=${atom_final})); then
      echo "bond     ${n}    $((${n}+1)) " >> ${InputPlaymolFile}
    fi
  done
  echo "" >> ${InputPlaymolFile}
done < structure.tmp
rm structure.tmp

cp ${xyz_file} xyzfile.temp
sed -e "1i $((${atom_count}-1))" -e \
       "1i # arquivo com as coordenadas xyz gerado a partir do script inp_play_homopdi" xyzfile.temp > ${xyz_file}
rm xyzfile.temp

echo "build    ${xyz_file}
box      density ${ro}

packmol  seed ${seed}   tolerance ${tol}   retry ${retry_value}  " >> ${InputPlaymolFile}
cat file1.tmp >> ${InputPlaymolFile}
rm file1.tmp
echo "packmol  action execute

write    lammps ${lammps_data_file} 
write    lammpstrj ${vmd_read_file}   " >> ${InputPlaymolFile}
echo -e  "**********Concluído**********\n"





# recuperando a forma do arquivo de entrada:
mv parameterfile.backup ${ParameterFile}

cd ../
done

echo "================SCRIPT PROCESSADO ATÉ O FINAL================"
