#!/bin/bash
#sort_into_categories $contig
#family member with the biggest number of not uniquely mapped reads reported

dir=$1;
contig=$2;

all_unmapped1=`grep $contig statistics/stat_mother | cut -f 3`; 
high_qual_unmapped1=`grep $contig statistics/stat_mother_high | cut -f 3`; 
diff_mother1=`expr $all_unmapped1 - $high_qual_unmapped1`;
percentage1=`echo "scale=3; $diff_mother1/$all_unmapped1*100" | bc`;

all_unmapped2=`grep $contig statistics/stat_father | cut -f 3`; 
high_qual_unmapped2=`grep $contig statistics/stat_father_high | cut -f 3`; 
diff_mother2=`expr $all_unmapped2 - $high_qual_unmapped2`;
percentage2=`echo "scale=3; $diff_mother2/$all_unmapped2*100" | bc`;

all_unmapped3=`grep $contig statistics/stat_daughter | cut -f 3`; 
high_qual_unmapped3=`grep $contig statistics/stat_daughter_high | cut -f 3`; 
diff_mother3=`expr $all_unmapped3 - $high_qual_unmapped3`;
percentage3=`echo "scale=3; $diff_mother3/$all_unmapped3*100" | bc`;

all_unmapped4=`grep $contig statistics/stat_son | cut -f 3`; 
high_qual_unmapped4=`grep $contig statistics/stat_son_high | cut -f 3`; 
diff_mother4=`expr $all_unmapped4 - $high_qual_unmapped4`;
percentage4=`echo "scale=3; $diff_mother4/$all_unmapped4*100" | bc`;

#searching for maximum - family member that has the biggest number of NOT uniquely mapped reads
    if (( $(echo "$percentage1 > $percentage2"|bc -l) ));
    	then
        	max1="$percentage1";
    	else
		max1="$percentage2";
    fi

    if (( $(echo "$percentage3 > $percentage4"|bc -l) ));
    	then
        	max2="$percentage3";
    	else
		max2="$percentage4";
    fi


    if (( $(echo "$max1 > $max2"|bc -l) ));
    	then
        	percentage="$max1";
    	else
		percentage="$max2";
    fi

echo $percentage;


