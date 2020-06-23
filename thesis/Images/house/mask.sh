for file in $(ls *mask.png); do
    python ../../../scripts/create_mask.py ../../../images/grey/bank.pgm $file
done
