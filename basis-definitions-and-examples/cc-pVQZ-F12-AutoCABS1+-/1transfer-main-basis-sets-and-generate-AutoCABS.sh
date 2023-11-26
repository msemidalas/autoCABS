cp cc-pVQZ-F12-orbital-basis-sets/*orca . 

for element in h he li be b c n o f ne na mg al si p s cl ar; do
    mkdir $element && mv $element.orca $element && cp AutoCABS1+-.py $element && cp myfunctions.py $element
done

#Generate basis set
for element in h he li be b c n o f ne na mg al si p s cl ar; do
    cd $element && python3 AutoCABS1+-.py && cd ..
done
cat /dev/null > cc-pVQZ-F12-AutoCABS1+-.mbas
for element in h he li be b c n o f ne na mg al si p s cl ar; do
    cat ${element}/${element}-AutoCABS.mpro >> cc-pVQZ-F12-AutoCABS1+-.mbas
done

sed -i '' 's/basis={//g;s/}//g' cc-pVQZ-F12-AutoCABS1+-.mbas

#SELF-DESTRUCT
#for element in h he li be b c n o f ne na mg al si p s cl ar;do rm -r $element;done 
