cp cc-pVTZ-F12-orbital-basis-sets/*orca . 

for element in h he li be b c n o f ne na mg al si p s cl ar; do
    mkdir $element && mv $element.orca $element && cp AutoCABS-generation-geom_mean_tight_and_diffuse_plus_single_layer.py $element && cp myfunctions.py $element
done

#Generate basis set
for element in h he li be b c n o f ne na mg al si p s cl ar; do
    cd $element && python3 AutoCABS-generation-geom_mean_tight_and_diffuse_plus_single_layer.py && cd ..
done
cat /dev/null > cc-pVTZ-F12-AutoCABS-geom_mean_tight_and_diffuse_plus_single_layer.mbas
for element in h he li be b c n o f ne na mg al si p s cl ar; do
    cat ${element}/${element}-AutoCABS.mpro >> cc-pVTZ-F12-AutoCABS-geom_mean_tight_and_diffuse_plus_single_layer.mbas
done

sed -i '' 's/basis={//g;s/}//g' cc-pVTZ-F12-AutoCABS-geom_mean_tight_and_diffuse_plus_single_layer.mbas

#SELF-DESTRUCT
#for element in h he li be b c n o f ne na mg al si p s cl ar;do rm -r $element;done 
