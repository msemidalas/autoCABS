cp cc-pVDZ-F12-orbital-basis-sets/*orca . 

for element in h he li be b c n o f ne na mg al si p s cl ar; do
    mkdir $element && mv $element.orca $element && cp AutoCABS2+-.py $element && cp myfunctions.py $element
done

#Generate basis set
for element in h he li be b c n o f ne na mg al si p s cl ar; do
    cd $element && python3 AutoCABS2+-.py && cd ..
done
cat /dev/null > cc-pVDZ-F12-AutoCABS2+-.mbas
for element in h he li be b c n o f ne na mg al si p s cl ar; do
    cat ${element}/${element}-AutoCABS.mpro >> cc-pVDZ-F12-AutoCABS2+-.mbas
done

sed -i '' 's/basis={//g;s/}//g' cc-pVDZ-F12-AutoCABS2+-.mbas

#SELF-DESTRUCT
#for element in h he li be b c n o f ne na mg al si p s cl ar;do rm -r $element;done 
