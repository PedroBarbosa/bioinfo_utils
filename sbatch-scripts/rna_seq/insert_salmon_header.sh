header="Name	Length	EffectiveLength	TPM	NumReads"

for file in "$@"
do
    sed -i -e "1i\
$header" $file
done
