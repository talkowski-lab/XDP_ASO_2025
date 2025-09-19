#!/bin/bash

# to get IR for all TAF1 exons
# run in data/
DIR=`pwd`

cd $DIR/IRFinder

for SAMPLE in `ls`
do 
	IR_1=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70586344"|grep "70587348"|cut -f20`
	echo -e "$SAMPLE\t$IR_1" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron1.txt
	IR_2=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70587463"|grep "70587903"|cut -f20`
	echo -e "$SAMPLE\t$IR_2" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron2.txt
	IR_3=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70588020"|grep "70595016"|cut -f20`
	echo -e "$SAMPLE\t$IR_3" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron3.txt
	IR_4=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70595136"|grep "70596799"|cut -f20`
	echo -e "$SAMPLE\t$IR_4" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron4.txt
	IR_5=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70597041"|grep "70597452"|cut -f20`
	echo -e "$SAMPLE\t$IR_5" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron5.txt
	IR_7=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70598303"|grep "70598673"|cut -f20`
	echo -e "$SAMPLE\t$IR_7" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron7.txt
	IR_8=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70598881"|grep "70601592"|cut -f20`
	echo -e "$SAMPLE\t$IR_8" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron8.txt
	IR_9=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70601769"|grep "70602385"|cut -f20`
	echo -e "$SAMPLE\t$IR_9" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron9.txt
	IR_10=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70602513"|grep "70602610"|cut -f20`
	echo -e "$SAMPLE\t$IR_10" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron10.txt
	IR_11=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70602718"|grep "70602840"|cut -f20`
	echo -e "$SAMPLE\t$IR_11" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron11.txt
	IR_12=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70603014"|grep "70603811"|cut -f20`
	echo -e "$SAMPLE\t$IR_12" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron12.txt
	IR_13=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70603985"|grep "70604794"|cut -f20`
	echo -e "$SAMPLE\t$IR_13" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron13.txt
	IR_14=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70604899"|grep "70607110"|cut -f20`
	echo -e "$SAMPLE\t$IR_14" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron14.txt
	IR_15=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70607311"|grep "70608086"|cut -f20`
	echo -e "$SAMPLE\t$IR_15" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron15.txt
	IR_17=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70608718"|grep "70609434"|cut -f20`
	echo -e "$SAMPLE\t$IR_17" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron17.txt
	IR_18=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70609515"|grep "70612418"|cut -f20`
	echo -e "$SAMPLE\t$IR_18" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron18.txt
	IR_19=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70612568"|grep "70612724"|cut -f20`
	echo -e "$SAMPLE\t$IR_19" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron19.txt
	IR_20=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70612844"|grep "70613150"|cut -f20`
	echo -e "$SAMPLE\t$IR_20" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron20.txt
	IR_21=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70613326"|grep "70613916"|cut -f20`
	echo -e "$SAMPLE\t$IR_21" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron21.txt
	IR_22=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70614095"|grep "70617102"|cut -f20`
	echo -e "$SAMPLE\t$IR_22" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron22.txt
	IR_23=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70617316"|grep "70618421"|cut -f20`
	echo -e "$SAMPLE\t$IR_23" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron23.txt
	IR_24=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70618587"|grep "70621377"|cut -f20`
	echo -e "$SAMPLE\t$IR_24" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron24.txt
	IR_25=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70621589"|grep "70626487"|cut -f20`
	echo -e "$SAMPLE\t$IR_25" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron25.txt
	IR_26=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70626596"|grep "70627423"|cut -f20`
	echo -e "$SAMPLE\t$IR_26" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron26.txt
	IR_27=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70627522"|grep "70627823"|cut -f20`
	echo -e "$SAMPLE\t$IR_27" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron27.txt
	IR_28=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70628001"|grep "70641158"|cut -f20`
	echo -e "$SAMPLE\t$IR_28" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron28.txt
	IR_29=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70641226"|grep "70642966"|cut -f20`
	echo -e "$SAMPLE\t$IR_29" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron29.txt
	IR_30=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70643089"|grep "70643823"|cut -f20`
	echo -e "$SAMPLE\t$IR_30" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron30.txt
	IR_31=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70643916"|grep "70644003"|cut -f20`
	echo -e "$SAMPLE\t$IR_31" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron31.txt
	IR_32a=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70644088"|grep "70664176"|cut -f20`
	echo -e "$SAMPLE\t$IR_32a" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron32a.txt
	IR_32all=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70644088"|grep "70676646"|cut -f20`
	echo -e "$SAMPLE\t$IR_32all" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron32all.txt
	IR_32b=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70664198"|grep "70676646"|cut -f20`
	echo -e "$SAMPLE\t$IR_32b" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron32b.txt
	IR_33=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70676714"|grep "70677217"|cut -f20`
	echo -e "$SAMPLE\t$IR_33" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron33.txt
	IR_34a=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70677334"|grep "70677483"|cut -f20`
	echo -e "$SAMPLE\t$IR_34a" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron34a.txt
	IR_34b=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70677334"|grep "70680717"|cut -f20`
	echo -e "$SAMPLE\t$IR_34b" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron34b.txt
	IR_35a=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70677489"|grep "70680717"|cut -f20`
	echo -e "$SAMPLE\t$IR_35a" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron35a.txt
	IR_35b=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70680834"|grep "70681626"|cut -f20`
	echo -e "$SAMPLE\t$IR_35b" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron35b.txt
	IR_36a=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70680843"|grep "70682028"|cut -f20`
	echo -e "$SAMPLE\t$IR_36a" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron36a.txt
	IR_36b=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70681728"|grep "70682028"|cut -f20`
	echo -e "$SAMPLE\t$IR_36b" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron36b.txt
	IR_37=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70682185"|grep "70683102"|cut -f20`
	echo -e "$SAMPLE\t$IR_37" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron37.txt
	IR_38a=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70683280"|grep "70686300"|cut -f20`
	echo -e "$SAMPLE\t$IR_38a" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron38a.txt
	IR_38b=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70683280"|grep "70707360"|cut -f20`
	echo -e "$SAMPLE\t$IR_38b" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron38b.txt
	IR_38c=`cat $SAMPLE/IRFinder-IR-dir.txt|grep "^X"|grep "70683280"|grep "70751018"|cut -f20`
	echo -e "$SAMPLE\t$IR_38c" >> $DIR/IRFinder_results/IR_TAF1/IR_Intron38c.txt
done
