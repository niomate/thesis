
for file in $(ls *inpaint*.pgm); do
    psnr=$(head -n25 $file | grep PSNR)
    echo $file
    echo $psnr
done
