cp cc-pVDZ-F12-orbital-basis-sets/*orca . 

for element in h he li be b c n o f ne na mg al si p s cl ar; do
    mkdir $element && mv $element.orca $element && cp AutoCABS-generation-geom_mean_tight_extended_2p_and_diffuse_plus_two_layers.py $element && cp myfunctions.py $element
done

#Generate basis set
for element in h he li be b c n o f ne na mg al si p s cl ar; do
    cd $element && python3 AutoCABS-generation-geom_mean_tight_extended_2p_and_diffuse_plus_two_layers.py && cd ..
done
cat /dev/null > cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
for element in h he li be b c n o f ne na mg al si p s cl ar; do
    cat ${element}/${element}-AutoCABS.mpro >> cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas
done

sed -i '' 's/basis={//g;s/}//g' cc-pVDZ-F12-AutoCABS-Geom_Mean_tight_extended_2p_and_diffuse_plus_two_layers.mbas

#SELF-DESTRUCT
#for element in h he li be b c n o f ne na mg al si p s cl ar;do rm -r $element;done 
