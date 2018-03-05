rm Errors.txt
rm Times.txt

echo "-1" | tr "\n" " " >> Errors.txt
echo "-1" | tr "\n" " " >> Times.txt


for j in {1..15}
do
  N=$[2**$j]
  echo $N | tr "\n" " " >> Errors.txt
  echo $N | tr "\n" " " >> Times.txt
done

for i in {1..5}
do
  Np=$[2**$i]
  echo "Np=$Np"
  echo "" >> Errors.txt
  echo "" >> Times.txt
  echo $Np | tr "\n" " " >> Errors.txt
  echo $Np | tr "\n" " " >> Times.txt

  for j in {1..15}
  do
    N=$[2**$j]
    echo "\t N=$N"
    make test N=$N Np=$Np >> tmp.txt
    tail -n 1 tmp.txt | tr "\n" " " >> Errors.txt
    { time make test N=$N Np=$Np ;} 2> tmp.txt
    head -n 2 tmp.txt | tail -n 1 |  tr -d "real ms" | tr -d " "| tr "\n" " " >> Times.txt
    rm tmp.txt
  done
  #echo "" >> Errors.txt
done
