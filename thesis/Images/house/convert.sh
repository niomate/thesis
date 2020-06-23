for file in $(ls *.pgm); do
    filename=$(basename -- "$file")
    extension="${filename##*.}"
    filename="${filename%.*}"
    convert $file "$filename.png"
done
