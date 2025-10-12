# Assembleur basé sur les graphes de Debruijn

Veuillez cloner le projet : 
```
git clone https://github.com/anaisdlss/tp_debruijn.git
```


Puis rendez vous dans le repertoire crée.
Créez un fichier *résultats* :
```
mkdir resultats
```

Puis exectutez le script :
```
python3 debruijn/debruijn.py -i data/eva71_plus_perfect.fq -k 22 -o resultats/contigs.fasta -f resultats/graph.png
```
Il est possible de modifier la taille des kmer avec -k


Après cela, créez la base de données BLAST pour la référence ```data/eva71.fna```:
```
makeblastdb -in data/eva71.fna -dbtype nucl
```
Puis comparez vos contigs obtenues dans ```resultats/contigs.fasta``` avec la référence :
```
blastn -query resultats/contigs.fasta -db data/eva71.fna -out resultats/blast_result.txt -outfmt 6
```

Les resultats obtenus se trouvent dans le fichier ```resultats/blast_result.txt```.
