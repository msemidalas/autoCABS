cp /gpfs01/work/emmanous/4F12basisforTMs/main-basis-set-definitions/cc-pVTZ-F12-orca/*orca .
for f in $(cat list-folders);do mkdir $f && mv $f.orca $f;done
for f in $(cat list-folders);do cp AutoCABS-generation-124.py $f;done
for f in $(cat list-folders);do cp myfunctions.py $f;done

#Generate basis set
for f in $(cat list-folders);do cd $f && python3 AutoCABS-generation-124.py && cd ..;done
cat /dev/null > cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat h/H-AutoCABS.mpro >> cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat he/He-AutoCABS.mpro >> cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat li/Li-AutoCABS.mpro >>  cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat be/Be-AutoCABS.mpro >>  cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat b/B-AutoCABS.mpro >>  cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat c/C-AutoCABS.mpro >>  cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat n/N-AutoCABS.mpro >>  cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat o/O-AutoCABS.mpro >>  cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat f/F-AutoCABS.mpro >>  cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat ne/Ne-AutoCABS.mpro >>  cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat na/Na-AutoCABS.mpro >> cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat mg/Mg-AutoCABS.mpro >> cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat al/Al-AutoCABS.mpro >> cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat si/Si-AutoCABS.mpro >> cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat p/P-AutoCABS.mpro >> cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat s/S-AutoCABS.mpro >> cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat cl/Cl-AutoCABS.mpro >> cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
cat ar/Ar-AutoCABS.mpro >> cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
sed -i 's/basis={//g' cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas
sed -i 's/}//g' cc-pVTZ-F12-AutoCABS-Geom_Mean_tight_and_diffuse_plus_single_layer.mbas

#SELF-DESTRUCT
#for f in $(cat list-folders);do rm -r $f;done 
