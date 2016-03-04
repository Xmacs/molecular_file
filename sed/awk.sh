for i in 20 50 100 200
do
	for j in 0.5 0 1 2.5 5 10 20
	do
		sed '1,10093d' ${i}_${j}.lammps > ${i}_${j}_1.lammps
		sed '100001,100022d' ${i}_${j}_1.lammps > ${i}_${j}_2.lammps
		awk '{ {print $2} }' ${i}_${j}_2.lammps > ./out/${i}_${j}.dat
		set a=0
		awk '{a+=$1;count+=1} END{print a/count}' ./out/${i}_${j}.dat > ./out/${i}_${j}_ave.dat
		rm -rf ${i}_${j}_1.lammps
		rm -rf ${i}_${j}_2.lammps
	done
done
