#collect Windows dependencies (DLLs) to include in the deployment
#dest_dir="./build/" #destination directory where to collect binary and DLLs
dest_dir=$1 #destination directory where to collect binary and DLLs
dlls_l=`objdump -x $dest_dir/wasa.exe |grep "DLL Name" | sed -e 's/^.* //g'` #find names of required DLLs
dlls_a=($dlls_l) #convert list to array
dlls_l=`echo "${dlls_a[@]}" | tr ' ' '\n' | sort | uniq | tr '\n' ' '` #remove duplicates
dlls_a=($dlls_l) #convert list to array



while :
do
    echo start loop 
	#echo DLL: "${dlls_a[@]}"
	# echo length: "${#dlls_a[@]}"
	echo DLL_l "$dlls_l"

	secondary_dependencies_l=$dlls_l #list of second-order dependencies, initialized with DLLs found so far
	for dll in $dlls_l
	do
	   echo ""
	   echo Dependencies of $dll
	   full_path=`which $dll`
	   echo ""
	   echo full path:$full_path
	   further_dlls=`objdump -x $full_path |grep "DLL Name" | sed -e 's/^.* //g'` #find names of required DLLs
	   echo ""
	   echo further: $further_dlls
	   #further_dlls=$further_dlls 
	   secondary_dependencies_l+=" $further_dlls" #collect further dependencies
	   #echo ""
	   #echo secondary_dependencies_l: $secondary_dependencies_l
	done   
	echo ""
	secondary_dependencies_a=($secondary_dependencies_l) #convert list to array
	secondary_dependencies_l=`echo "${secondary_dependencies_a[@]}" | xargs | tr ' ' '\n'  | sort | uniq | tr '\n' ' '` #remove duplicates
	secondary_dependencies_a=($secondary_dependencies_l) #convert list to array
		
	# echo DLL: "${dlls_a[@]}"
	# echo length: "${#dlls_a[@]}"
	# echo DLL_l "$dlls_l"
	# echo "$dlls_l" > dll.txt
#	read -n 1 KEY
	# echo secondary_dependencies_a: "${secondary_dependencies_a[@]}"
	# echo length: "${#secondary_dependencies_a[@]}"
	# echo secondary_dependencies_l: "$secondary_dependencies_l"
	# echo "$secondary_dependencies_l" > sec.txt
	# read -n 1 KEY

	if [ "$dlls_l" == "$secondary_dependencies_l" ]
	then
		echo same
		break
	else
		echo different, next loop
	    dlls_l=$secondary_dependencies_l 
		dlls_a=($dlls_l) #convert list to array
	fi
done

echo DLLs to copy: "$dlls_l"

for dll in $dlls_l
do
   echo $dll
   full_path=`which $dll`
   echo $full_path
   if echo $full_path | grep Windows ; then
     echo "is win"
   else
     echo "copy"
     cp $full_path $dest_dir #not in windows directory, copy to collect
   fi
done
