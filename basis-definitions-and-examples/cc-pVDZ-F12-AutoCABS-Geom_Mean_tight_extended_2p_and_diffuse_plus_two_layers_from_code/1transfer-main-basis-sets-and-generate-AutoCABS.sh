cp /gpfs01/work/emmanous/4F12basisforTMs/main-basis-set-definitions/cc-pVDZ-F12-orca/*orca .
for f in $(cat list-folders);do mkdir $f && mv $f.orca $f;done
for f in $(cat list-folders);do cp AutoCABS-generation-225-tight-extended-2p.py $f;done
for f in $(cat list-folders);do cp myfunctions.py $f;done

#Generate basis set
for f in $(cat list-folders);do cd $f && python3 AutoCABS-generation-225-tight-extended-2p.py && cd ..;done
cat /dev/null > cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat h/H-AutoCABS.mpro >> cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat he/He-AutoCABS.mpro >> cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat li/Li-AutoCABS.mpro >>  cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat be/Be-AutoCABS.mpro >>  cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat b/B-AutoCABS.mpro >>  cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat c/C-AutoCABS.mpro >>  cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat n/N-AutoCABS.mpro >>  cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat o/O-AutoCABS.mpro >>  cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat f/F-AutoCABS.mpro >>  cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat ne/Ne-AutoCABS.mpro >>  cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat na/Na-AutoCABS.mpro >> cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat mg/Mg-AutoCABS.mpro >> cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat al/Al-AutoCABS.mpro >> cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat si/Si-AutoCABS.mpro >> cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat p/P-AutoCABS.mpro >> cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat s/S-AutoCABS.mpro >> cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat cl/Cl-AutoCABS.mpro >> cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
cat ar/Ar-AutoCABS.mpro >> cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
sed -i 's/basis={//g' cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
sed -i 's/}//g' cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas

#SELF-DESTRUCT
#for f in $(cat list-folders);do rm -r $f;done 
