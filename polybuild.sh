#!/bin/bash
#=========================================================================================================
# polybuild
#      Program that receives an input file with information from a dispersed linear 
#      architecture polymeric microstructure and writes an input file to be processed 
#      by Playmol to package the microstructure generating an initial configuration
#      for use in molecular dynamics simulations.
#
# Developed by Tiago S. M. Lemos
#      Programa de Engenharia Química/COPPE
#      Universidade Federal do Rio de Janeiro – CP: 69502
#      Cidade Universitária, Rio de Janeiro, RJ, 21941-972, Brazil
#      e-mail: tsmlemos@hotmail.com.br
# 
# Rio de Janeiro - October 1, 2016
#========================================================================================================
# If this script is used in the preparation of scientific publications, the following paper should be cited:
#
# LEMOS, T. S. M., ABREU, C. R. A., PINTO, J. C. 
# Mesoscopic Simulation of Dispersed Copolymers - Effects of Chain Length, Chemical
#    Composition, and Block Length Distributions on Self-Assembly
# Macromolecular Theory and Simulations, 2019
#========================================================================================================
# Recognize the input file:

for ParameterFile in $*
do
echo "polybuild is processing the \"${ParameterFile}\" input file"
cp ${ParameterFile} parameterfile.backup
sed 's/[ ]\+/ /g ; /^ *\?$/d' parameterfile.backup > ${ParameterFile}
prefix=$(echo "${ParameterFile}" | cut -d. -f 1)


#========================================================================================================
# Define the output files names:

InputPlaymolFile="${prefix}.mol"                  # Playmol input file name
xyz_file="${prefix}.xyz"                          # xyz file name
lammpstrj_read_file="${prefix}.lammpstrj"         # Lammps trajectory file name
lammps_data_file="data.${prefix}.lmp"             # name of data file to generate for Lammps simulation


#========================================================================================================
# Store the parameter values for the system to be simulated:

k_spring=$(            grep -E "^k_spring "               ${ParameterFile}  | cut -d " " -f 3)
xeq_spring=$(          grep -E "^xeq_spring "             ${ParameterFile}  | cut -d " " -f 3)


#========================================================================================================
# Store Playmol packing parameters:

ro=$(                  grep -E "^ro "                   ${ParameterFile}  | cut -d " " -f 3)
r_bond=$(              grep -E "^r_bond "               ${ParameterFile}  | cut -d " " -f 3)
seed=$(                grep -E "^seed "                 ${ParameterFile}  | cut -d " " -f 3)
tol=$(                 grep -E "^tol "                  ${ParameterFile}  | cut -d " " -f 3)
retry_value=$(         grep -E "^retry_value "          ${ParameterFile}  | cut -d " " -f 3)
                                                                                                               

#========================================================================================================
# Recognize the types of atoms/beads and make a list with their names:

init_micro=$(($(    grep -E -n "Microstructure"                  ${ParameterFile}  | cut -d ":" -f 1)+1))
final_micro=$(($(   grep -E -n "List of components"              ${ParameterFile}  | cut -d ":" -f 1)-1)) 
init_compmass=$(($( grep -E -n "List of components"              ${ParameterFile}  | cut -d ":" -f 1)+1)) 
final_compmass=$(($(grep -E -n "Definition of system parameters" ${ParameterFile}  | cut -d ":" -f 1)-1)) 

sed -n ${init_micro},${final_micro}p        ${ParameterFile} > structure.tmp
sed -n ${init_compmass},${final_compmass}p  ${ParameterFile} > components.tmp


#========================================================================================================
# Write the Playmol input file:

echo -e "**********Writing input for Playmol********** \n. \n.. \n... "
echo "# File generated from 'polybuild' script
# creation date: $(date) 

  "  > ${InputPlaymolFile}

while read id type mass; do
  echo "atom_type       ${type} 
mass            ${type}       ${mass}"  >> ${InputPlaymolFile}
done < components.tmp
rm components.tmp

echo "define          ro  as ${ro}            # numeric density
bond_type       *   *   ${k_spring} ${xeq_spring}
  
"  >> ${InputPlaymolFile}


# Read the microstructure and generate a configuration file in xyz format:

atom_count=1
molecule_count=1
> file1.tmp
> ${xyz_file}
while read line; do

  echo "# molecule  ${molecule_count} " >> ${InputPlaymolFile}
  n_words="$(echo "${line}" | wc -w)"

  quantity="$(echo "${line}" | cut -d " " -f ${n_words})"
  echo "# N${molecule_count} =  ${quantity}" >> ${InputPlaymolFile}
  echo "packmol    copy ${molecule_count}  ${quantity} " >> file1.tmp
  molecule_count=$((${molecule_count}+1))

  atom_init=${atom_count}
  x_atom=0.000 

  for n in $(seq 1 2 $((${n_words}-1))) ; do
    type="$(echo "${line}" | cut -d " " -f ${n})"
    block_size="$(echo "${line}" | cut -d " " -f $((${n}+1)))"

    for n2 in $(seq 1 ${block_size}); do 
      echo "atom     ${atom_count}    ${type}" >> ${InputPlaymolFile}
      echo "${atom_count}      ${x_atom}    0.000    0.000" >> ${xyz_file}
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
       "1i # file with xyz coordinates generated from polybuild script" xyzfile.temp > ${xyz_file}
rm xyzfile.temp

echo "build    ${xyz_file}
box      density ${ro}

packmol  seed ${seed}   tolerance ${tol}   retry ${retry_value}  " >> ${InputPlaymolFile}
cat file1.tmp >> ${InputPlaymolFile}
rm file1.tmp
echo "packmol  action execute

write    lammps ${lammps_data_file} 
write    lammpstrj ${lammpstrj_read_file}   " >> ${InputPlaymolFile}
echo -e  "**********Completed**********\n"


#========================================================================================================
# Restore input file format:
mv parameterfile.backup ${ParameterFile}


#========================================================================================================
done


