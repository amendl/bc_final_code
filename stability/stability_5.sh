cd ~
cd work/2/NH3_psi4_HF
python /home/users/mendla/work/2/psi4_analysis/psi4_higher_stability.py --files="AVDZ/dma4.punch;VDZ/dma4.punch;AVTZ/dma4.punch;VTZ/dma4.punch;AVQZ/dma4.punch;VQZ/dma4.punch;AV5Z/dma4.punch;V5Z/dma4.punch" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --site_index="0" --output="NH3_high_H"
cd ~
cd work/2/NH3_psi4
python /home/users/mendla/work/2/psi4_analysis/psi4_higher_stability.py --files="AVDZ/dma4.punch;VDZ/dma4.punch;AVTZ/dma4.punch;VTZ/dma4.punch;AVQZ/dma4.punch;VQZ/dma4.punch;AV5Z/dma4.punch;V5Z/dma4.punch" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --site_index="0" --output="NH3_high_CC"
cd ~
cd work/2/HF_psi4_HF
python /home/users/mendla/work/2/psi4_analysis/psi4_higher_stability.py --files="AVDZ/dma4.punch;VDZ/dma4.punch;AVTZ/dma4.punch;VTZ/dma4.punch;AVQZ/dma4.punch;VQZ/dma4.punch;AV5Z/dma4.punch;V5Z/dma4.punch" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --site_index="1" --output="HF_high_H"
cd ~
cd work/2/HF_psi4
python /home/users/mendla/work/2/psi4_analysis/psi4_higher_stability.py --files="AVDZ/dma4.punch;VDZ/dma4.punch;AVTZ/dma4.punch;VTZ/dma4.punch;AVQZ/dma4.punch;VQZ/dma4.punch;AV5Z/dma4.punch;V5Z/dma4.punch" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --site_index="1" --output="HF_high_CC"
cd ~
cd work/2/CO2_psi4_HF
python /home/users/mendla/work/2/psi4_analysis/psi4_higher_stability.py --files="AVDZ/dma4.punch;VDZ/dma4.punch;AVTZ/dma4.punch;VTZ/dma4.punch;AVQZ/dma4.punch;VQZ/dma4.punch;AV5Z/dma4.punch;V5Z/dma4.punch" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --site_index="1" --output="CO2_high_H"
cd ~
cd work/2/CO2_psi4
python /home/users/mendla/work/2/psi4_analysis/psi4_higher_stability.py --files="AVDZ/dma4.punch;VDZ/dma4.punch;AVTZ/dma4.punch;VTZ/dma4.punch;AVQZ/dma4.punch;VQZ/dma4.punch;AV5Z/dma4.punch;V5Z/dma4.punch" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --site_index="1" --output="CO2_high_CC"
cd ~
cd work/2/methanol_psi4_HF
python /home/users/mendla/work/2/psi4_analysis/psi4_higher_stability.py --files="AVDZ/dma4.punch;VDZ/dma4.punch;AVTZ/dma4.punch;VTZ/dma4.punch;AVQZ/dma4.punch;VQZ/dma4.punch;AV5Z/dma4.punch;V5Z/dma4.punch" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ;AV5Z;V5Z" --site_index="0" --output="methanol_high_H"
cd ~
cd work/2/methanol_psi4
python /home/users/mendla/work/2/psi4_analysis/psi4_higher_stability.py --files="AVDZ/dma4.punch;VDZ/dma4.punch;AVTZ/dma4.punch;VTZ/dma4.punch;AVQZ/dma4.punch;VQZ/dma4.punch" --labels="AVDZ;VDZ;AVTZ;VTZ;AVQZ;VQZ" --site_index="0"  --output="methanol_high_CC"


