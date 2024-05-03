# cd ~
# cd work/2/HF2
# python /home/users/mendla/work/stability.py --files="01_AVDZ/molpro.out;01_VDZ/molpro.out;01_AVTZ/molpro.out;01_VTZ/molpro.out;01_AVQZ/molpro.out;01_VQZ/molpro.out;01_AV5Z/molpro.out;01_V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="3" --site_index="0" --no_of_records="5" --output="HF_H"
cd ~
cd work/1/13
python /home/users/mendla/work/stability.py --files="AVDZ/molpro.out;VDZ/molpro.out;AVTZ/molpro.out;VTZ/molpro.out;AVQZ/molpro.out;VQZ/molpro.out;AV5Z_simple/molpro.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="0" --no_of_records="1" --output="methanol_H"
# cd ~
python /home/users/mendla/work/stability.py --files="AVDZ/molpro_cc.out;VDZ/molpro_cc.out;AVTZ/molpro_cc.out;VTZ/molpro_cc.out;AVQZ/molpro_cc.out;VQZ/molpro_cc.out;AV5Z_simple/molpro_cc.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="0" --no_of_records="1" --output="methanol_CC"
cd ~
cd work/2/CO2_molpro_HF
python /home/users/mendla/work/stability.py --files="AVDZ/molpro.out;VDZ/molpro.out;AVTZ/molpro.out;VTZ/molpro.out;AVQZ/molpro.out;VQZ/molpro.out;AV5Z/molpro.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="1" --no_of_records="1" --output="CO2_H"
cd ~
cd work/2/CO2_molpro
python /home/users/mendla/work/stability.py --files="AVDZ/molpro.out;VDZ/molpro.out;AVTZ/molpro.out;VTZ/molpro.out;AVQZ/molpro.out;VQZ/molpro.out;AV5Z/molpro.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="1" --no_of_records="1" --output="CO2_CC"
cd ~
cd work/2/HF_molpro_HF
python /home/users/mendla/work/stability.py --files="AVDZ/molpro.out;VDZ/molpro.out;AVTZ/molpro.out;VTZ/molpro.out;AVQZ/molpro.out;VQZ/molpro.out;AV5Z/molpro.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="1" --no_of_records="1" --output="HF_H"
cd ~
cd work/2/HF_molpro
python /home/users/mendla/work/stability.py --files="AVDZ/molpro.out;VDZ/molpro.out;AVTZ/molpro.out;VTZ/molpro.out;AVQZ/molpro.out;VQZ/molpro.out;AV5Z/molpro.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="1" --no_of_records="1" --output="HF_CC"
cd ~
cd work/2/NH3_molpro_HF
python /home/users/mendla/work/stability.py --files="AVDZ/molpro.out;VDZ/molpro.out;AVTZ/molpro.out;VTZ/molpro.out;AVQZ/molpro.out;VQZ/molpro.out;AV5Z/molpro.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="0" --no_of_records="1" --output="NH3_H"
cd ~
cd work/2/NH3_molpro
python /home/users/mendla/work/stability.py --files="AVDZ/molpro.out;VDZ/molpro.out;AVTZ/molpro.out;VTZ/molpro.out;AVQZ/molpro.out;VQZ/molpro.out;AV5Z/molpro.out;V5Z/molpro.out" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --record_index="0" --site_index="0" --no_of_records="1" --output="NH3_CC"
