for e in 1000 ; do 
for B in 5 10 20 50 ; do  
for n in 5 10 20 50 ; do  
    python -u dpboot_composition.py -n $n -B $B -e $e &
done
done
done;