# Spawns a process for each index to compress using WAH

for f in `ls indexes/`
do
    ~/Python-3.2.3/python RLEindex.py < indexes/$f > $f.wahsize &
done

for job in `jobs -p`
do
    wait $job
done

