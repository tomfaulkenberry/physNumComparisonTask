subject_nr=$(cat ~/Desktop/mousetrackingOperatorPreviewEffect/procedure/subject_nr.txt)
next_subject=$(($subject_nr + 1))

echo "Subject number = $subject_nr"
echo $next_subject > ~/Desktop/mousetrackingOperatorPreviewEffect/procedure/subject_nr.txt

if [ ! -d data ]; then
	mkdir data
fi


if ((subject_nr % 2 == 0)); then
	echo "Starting experiment..."
	opensesamerun ~/Desktop/mousetrackingOperatorPreviewEffect/procedure/Experiment-TF.osexp -s $subject_nr -l "/home/mathcog/Desktop/mousetrackingOperatorPreviewEffect/results/data/subject_$subject_nr.csv" -f

else
	echo "Starting experiment..."
	opensesamerun ~/Desktop/mousetrackingOperatorPreviewEffect/procedure/Experiment-FT.osexp -s $subject_nr -l "/home/mathcog/Desktop/mousetrackingOperatorPreviewEffect/results/data/subject_$subject_nr.csv" -f

fi



echo Experiment complete!
echo Goodbye.
