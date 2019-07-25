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

echo "processando o integ_playlammps para o sistema ${ParameterFile}"
#ParameterFile=$1


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


# Leia as flags do arquivo de configuração para as instruções:

Playmol=$(       grep -E "^Playmol "            ${ParameterFile}  | cut -d " " -f 3)
Manter_xyzfile=$(grep -E "^Manter_xyzfile "     ${ParameterFile}  | cut -d " " -f 3)
Lammps=$(        grep -E "^Lammps "             ${ParameterFile}  | cut -d " " -f 3)
Sq_file=$(       grep -E "^Sq_file "             ${ParameterFile}  | cut -d " " -f 3)


# Defina os nomes dos arquivos de saída:

if [ "${Playmol}" = "yes" ]; then
  InputPlaymolFile="${prefixo}.mol"                  # nome do arquivo de input a ser lido pelo Playmol
  xyz_file="${prefixo}.xyz"                          # nome do arquivo xyz a ser gerado
  vmd_read_file="${prefixo}.lammpstrj"               # nome do arquivo passível de visualização por VMD
fi

if [ "${Lammps}" = "yes" ]; then
  lammps_input="${prefixo}.lmp"                      # nome do arquivo de input da simulação Lammps
  log_file="${prefixo}.log"                          # nome do arquivo de output da simulação Lammps
fi
lammps_data_file="data.${prefixo}.lmp"               # nome do arquivo de dados gerado para a simulação Lammps 


# Definição das variáveis do sistema a ser simulado:

k_spring=$(            grep -E "^k_spring "               ${ParameterFile}  | cut -d " " -f 3)
xeq_spring=$(          grep -E "^xeq_spring "             ${ParameterFile}  | cut -d " " -f 3)


if [ "${Playmol}" = "yes" ]; then
  # Definição dos parâmetros Playmol:                                                                  
  
  ro=$(                  grep -E "^ro "                   ${ParameterFile}  | cut -d " " -f 3)
  r_bond=$(              grep -E "^r_bond "               ${ParameterFile}  | cut -d " " -f 3)
                                                                                                       
  seed=$(                grep -E "^seed "                 ${ParameterFile}  | cut -d " " -f 3)
  tol=$(                 grep -E "^tol "                  ${ParameterFile}  | cut -d " " -f 3)
  retry_value=$(         grep -E "^retry_value "          ${ParameterFile}  | cut -d " " -f 3)
fi                                                      
                                                                                                               

if [ "${Lammps}" = "yes" ]; then
  # Parâmetros para simulação DPD (LAMMPS):                                                                      

  dimension=$(           grep -E "^dimension "            ${ParameterFile}  | cut -d " " -f 3)
  atom_style=$(          grep -E "^atom_style "           ${ParameterFile}  | cut -d " " -f 3)
  units=$(               grep -E "^units "                ${ParameterFile}  | cut -d " " -f 3)
  dt=$(                  grep -E "^dt "                   ${ParameterFile}  | cut -d " " -f 3)
  steps=$(               grep -E "^steps "                ${ParameterFile}  | cut -d " " -f 3)
  steps_final=$(         grep -E "^steps_final "          ${ParameterFile}  | cut -d " " -f 3)
  qui=$(                 grep -E "^qui "                  ${ParameterFile}  | cut -d " " -f 3)
  bondstyle=$(           grep -E "^bondstyle "            ${ParameterFile}  | cut -d " " -f 3)                             
  Temp=$(                grep -E "^Temp "                 ${ParameterFile}  | cut -d " " -f 3)
  r_cut=$(               grep -E "^r_cut "                ${ParameterFile}  | cut -d " " -f 3)
  gamma=$(               grep -E "^gamma "                ${ParameterFile}  | cut -d " " -f 3)
  steps_random=$(        grep -E "^steps_random "         ${ParameterFile}  | cut -d " " -f 3)
  steps_restart=$(       grep -E "^steps_restart "        ${ParameterFile}  | cut -d " " -f 3)
  seed1=$(               grep -E "^seed1 "                ${ParameterFile}  | cut -d " " -f 3)
  seed2=$(               grep -E "^seed2 "                ${ParameterFile}  | cut -d " " -f 3)
  special_bonds_x=$(     grep -E "^special_bonds_x "      ${ParameterFile}  | cut -d " " -f 3)
  special_bonds_y=$(     grep -E "^special_bonds_y "      ${ParameterFile}  | cut -d " " -f 3)
  special_bonds_z=$(     grep -E "^special_bonds_z "      ${ParameterFile}  | cut -d " " -f 3)
  comm_modify_cutoff=$(  grep -E "^comm_modify_cutoff "   ${ParameterFile}  | cut -d " " -f 3)
  neighbor_bin=$(        grep -E "^neighbor_bin "         ${ParameterFile}  | cut -d " " -f 3)
  neigh_modify_every=$(  grep -E "^neigh_modify_every "   ${ParameterFile}  | cut -d " " -f 3)
  neigh_modify_delay=$(  grep -E "^neigh_modify_delay "   ${ParameterFile}  | cut -d " " -f 3)
  neigh_modify_cluster=$(grep -E "^neigh_modify_cluster " ${ParameterFile}  | cut -d " " -f 3)
  freq_ppm=$(            grep -E "^freq_ppm "             ${ParameterFile}  | cut -d " " -f 3)
  freq_lammpstrj=$(      grep -E "^freq_lammpstrj "       ${ParameterFile}  | cut -d " " -f 3)
  freq_thermo=$(         grep -E "^freq_thermo "          ${ParameterFile}  | cut -d " " -f 3)

  # Leitura dos potenciais de pares.
  init_potpair=$(($( grep -E -n "Coeficiente de interação" ${ParameterFile}  | cut -d ":" -f 1)+1)) 
  final_potpair=$(($(grep -E -n "Definição das variáveis"  ${ParameterFile}  | cut -d ":" -f 1)-1)) 

  sed -n ${init_potpair},${final_potpair}p       ${ParameterFile} > potpairs.tmp
fi

if [ "${Playmol}" = "yes" ]; then
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
  
  
  # Invocando o Playmol para a leitura do script gerado:
  
  echo -e "**********Invocando o Playmol********** \n. \n.. \n... "
  playmol -i ${InputPlaymolFile} > out.playmol
  rm out.playmol
  if [ "${Manter_xyzfile}" = "no" ]; then
    rm  ${xyz_file}
  fi
  echo -e  "**********Concluído**********\n"
  
 
  # formatar o arquivo 'data' gerado:
  
  cp ${lammps_data_file} data.tmp
  #########################################################################################
  # (ATUALIZADO: 07/02/19)
  #numero_comp=$((${final_compmass}-${init_compmass}+1))
  #sed '/Pair Coeffs/,/Atoms/d;/Angles/,$d' data.tmp | \
  #      sed -e "$(($( sed '/Masses/q' data.tmp   | grep -En Masses | cut -d: -f1)+3+${numero_comp}))i Atoms" > ${lammps_data_file}
  sed '/Atoms/p'  data.tmp | sed '/\# Elements/,/Atoms/ d' > ${lammps_data_file}
  #########################################################################################
  rm data.tmp

fi


if [ "${Lammps}" = "yes" ]; then
  # Escrita do input Lammps:
  
  echo "log ${log_file}

# Arquivo de input para a simulação de ${prefixo}

dimension       ${dimension}
atom_style      ${atom_style}
units           ${units}

read_data       ${lammps_data_file}
# read_restart  restart.file1

# Lista de parâmetros da simulação

variable seed1                equal   ${seed1}
variable seed2                equal   ${seed2}

variable dt                   equal   ${dt}        # passo temporal
variable steps                equal   ${steps}     # quantidade de passos de tempos 
variable steps_final          equal   ${steps_final}   # quantidade de passos de tempos para cálculo de S(q)

variable qui                  equal   ${qui}  " > ${lammps_input} 

  while read i j a_ij; do
    echo "variable a_${i}${j}                 equal   ${a_ij} " >> ${lammps_input}
  done < potpairs.tmp
 
  echo "
variable bondstyle            equal   ${bondstyle}
variable k_spring             equal   ${k_spring}
variable xeq_spring           equal   ${xeq_spring}

variable Temp                 equal   ${Temp}
variable r_cut                equal   ${r_cut}
variable gamma                equal   ${gamma}

variable special_bonds_x      equal   ${special_bonds_x}
variable special_bonds_y      equal   ${special_bonds_y}
variable special_bonds_z      equal   ${special_bonds_z}
                                     
variable comm_modify_cutoff   equal   ${comm_modify_cutoff}              
variable neighbor_bin         equal   ${neighbor_bin}          
variable neigh_modify_every   equal   ${neigh_modify_every}
variable neigh_modify_delay   equal   ${neigh_modify_delay}
variable neigh_modify_cluster equal   ${neigh_modify_cluster}
                                     
variable steps_random         equal   ${steps_random}
variable steps_restart        equal   ${steps_restart}
                                     
variable freq_lammpstrj       equal   ${freq_lammpstrj}
variable freq_ppm             equal   ${freq_ppm}
variable freq_thermo          equal   ${freq_thermo}


#=========================================================================

#group         typeA  type 1
#group         typeB  type 2

#bond_style    \${bondstyle}
bond_style    ${bondstyle}
bond_coeff    * \${k_spring} \${xeq_spring}

# pair_style dpd  T  cutoff seed
pair_style dpd  \${Temp} \${r_cut}  \${seed1}

# pair_coeff  I  J     A      gamma     cutoff " >> ${lammps_input} 

  while read i j a_ij; do
    echo "pair_coeff  ${i}  ${j}    \${a_${i}${j}}  \${gamma}  \${r_cut} "  >> ${lammps_input}
  done < potpairs.tmp
  

  echo "
special_bonds lj \${special_bonds_x} \${special_bonds_y} \${special_bonds_z}
comm_modify   vel yes cutoff \${comm_modify_cutoff}
#neighbor      \${neighbor_bin} bin
neigh_modify  every \${neigh_modify_every} delay \${neigh_modify_delay} cluster ${neigh_modify_cluster}

velocity      all  create \${Temp} \${seed2} dist gaussian
timestep      \${dt}
fix           1  all nve  

#======================================================================== 
# randomização da geometria de partida:" >> ${lammps_input}

  while read i j a_ij; do
    echo "pair_coeff  ${i}  ${j}    \${a_11}  \${gamma}  \${r_cut} "  >> ${lammps_input}
  done < potpairs.tmp

  #pair_coeff  1  2  \${a_ii}  \${gamma}  \${r_cut}
  echo "thermo \${freq_thermo} " >> ${lammps_input}
  echo "run \${steps_random} " >> ${lammps_input}
  echo "reset_timestep 0 " >> ${lammps_input}
  
  while read i j a_ij; do
    echo "pair_coeff  ${i}  ${j}    \${a_${i}${j}}  \${gamma}  \${r_cut} "  >> ${lammps_input}
  done < potpairs.tmp
  rm potpairs.tmp
  
  #pair_coeff  1  2  \${a_ij}  \${gamma}  \${r_cut}
  
  
  echo " 
#======================================================================== 
# simulação:

# velocity group-ID style args keyword value ... 
dump 2 all image \${freq_ppm} dump_${prefixo}.*.ppm type type
run \${steps}

dump 1 all atom \${freq_lammpstrj} trjsimul_${prefixo}.lammpstrj 
run \${steps_final}
clear" >> ${lammps_input}
fi


if [ "${Sq_file}" = "yes" ]; then
  # Leitura dos Parâmetros

  inputfileformat=$(         grep -E "^inputfile_format "           ${ParameterFile}  | cut -d " " -f 3)
  partial_Sq=$(              grep -E "^partial_Sq "                 ${ParameterFile}  | cut -d " " -f 3)
  number_of_comps_to_Sq=$(   grep -E "^number_of_components_to_Sq " ${ParameterFile}  | cut -d " " -f 3)
  name_of_comps_to_Sq=$(     grep -E "^name_of_components "         ${ParameterFile}  | cut -d "=" -f 2 | cut -d "#" -f 1 )
  equilibrations=$(          grep -E "^number_of_equilibrations "   ${ParameterFile}  | cut -d " " -f 3)
  configurations=$(          grep -E "^number_of_configurations_to_sample "   ${ParameterFile}  | cut -d " " -f 3)
  wave_vector_x=$(           grep -E "^wave_vector_x_dimension "    ${ParameterFile}  | cut -d " " -f 3)
  wave_vector_y=$(           grep -E "^wave_vector_y_dimension "    ${ParameterFile}  | cut -d " " -f 3)
  wave_vector_z=$(           grep -E "^wave_vector_z_dimension "    ${ParameterFile}  | cut -d " " -f 3)
  tolerance=$(               grep -E "^tolerance "                  ${ParameterFile}  | cut -d " " -f 3)
  
  # Escrita do script para o cômputo de S(q) Lammps:
 
  echo "trjsimul_${prefixo}.lammpstrj ! name of input file
${inputfileformat} ! input file format (0=xyz ; 1=lammpstrj)
${prefixo} ! name of system
                                      
.${partial_Sq}. ! calculate S(q) global (.false.) or partial (.true.)
${number_of_comps_to_Sq} ! numbers of components to S(q) calculation
${name_of_comps_to_Sq} ! name of components to partial structure factor calculation
                      
${equilibrations} ! number of equilibrations
${configurations} ! number of configurations to be sampled
${wave_vector_x} ! number of wave vector for x dimension
${wave_vector_y} ! number of wave vector for y dimension
${wave_vector_z} ! number of wave vector for z dimension
${tolerance} ! tolerance on sampler os S(q) and q
                                      
if (infile_format == 0) then                   
  -6.94                    ! x minimum value on xyz file
  6.94                     ! \" maximum   \"   \"   \"   \"
  -6.94                    ! y minimum   \"   \"   \"   \"
  6.94                     ! \" maximum   \"   \"   \"   \"
  -6.94                    ! z minimum   \"   \"   \"   \"
  6.94                     ! \" maximum   \"   \"   \"   \" 
end if" > sq_${prefixo}
fi


# recuperando a forma do arquivo de entrada:

#mv parameterfile.backup $1
mv parameterfile.backup ${ParameterFile}

cd ../
done


### Temporário: CORREÇÃO DO PLAYMOL
# (ATUALIZAÇÃO: PROBLEMA CORRIGIDO EM NOVA VERSÃO - 07/01/2019) 
#echo "Complementação temporária do software:"
#corrige_playmol $1
#corrige_playmol ${ParameterFile}

echo "================SCRIPT PROCESSADO ATÉ O FINAL================"
