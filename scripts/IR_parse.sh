#/bin/bash

# scripts after running IRfinder, adapted from Rachita
# run in data

DIR=`pwd`

cd $DIR/IRFinder

for SAMPLE in `ls`
do 
  IR=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70644088"|grep "70676646"|cut -f20`
  echo -e "$SAMPLE\t$IR" >> $DIR/IRFinder_results/IR_Intron32.txt
done

for SAMPLE in `ls`
do 
  SplicedReadsProximalToSVA=`cat $SAMPLE/IRFinder-JuncCount.txt|grep "^X"|grep "70644088"|grep "70659646"|cut -f7`
  Ex32_Ex32splice=`cat $SAMPLE/IRFinder-JuncCount.txt|grep "^X"|grep "70644088"|grep "70676646"|cut -f7`
  echo -e "$SAMPLE\t$SplicedReadsProximalToSVA\t$Ex32_Ex32splice" >> $DIR/IRFinder_results/i32_junctionSpliced.txt
done

for SAMPLE in `ls`
do 
  IR=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70644088"|grep "70659646"|cut -f19`
  Ex32_33=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70644088"|grep "70676646"|cut -f19`
  echo -e "$SAMPLE\t$IR\t$Ex32_33" >> $DIR/IRFinder_results/AS_inIntron32.txt
done

for SAMPLE in `ls`
do 
  SplicedReadsProximalToSVA=`cat $SAMPLE/IRFinder-JuncCount.txt|grep "^X"|grep "70644088"|grep "70659646" | cut -f7`
  Ex32_Ex32splice=`cat $SAMPLE/IRFinder-JuncCount.txt|grep "^X"|grep "70644088"|grep "70676646"|cut -f7`
  echo -e "$SAMPLE\t$SplicedReadsProximalToSVA\t$Ex32_Ex32splice" >> $DIR/IRFinder_results/SpliceqPCRData.txt
done

for SAMPLE in `ls`
do 
  SplicedReadsProximalToSVA=`cat $SAMPLE/IRFinder-JuncCount.txt|grep "^X"|grep "70644088"|awk '{if($3 < 70660363){sum+=$5}{ print sum }}' |  tail -1`
  Ex32_Ex32splice=`cat $SAMPLE/IRFinder-JuncCount.txt|grep "^X"|grep "70644088"|grep "70676646"|cut -f7`
  echo -e "$SAMPLE\t$SplicedReadsProximalToSVA\t$Ex32_Ex32splice" >> $DIR/IRFinder_results/AltSpliced_proximalToSVA.txt
done

