#collect Windows dependencies (DLLs) to include in the deployment
#dest_dir="./build/" #destination directory where to collect binary and DLLs
dest_dir=$1 #destination directory where to collect binary and DLLs
dlls=`objdump -x a.exe |grep "DLL Name" | sed -e 's/^.* //g'` #find names of required DLLs
for dll in $dlls
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
