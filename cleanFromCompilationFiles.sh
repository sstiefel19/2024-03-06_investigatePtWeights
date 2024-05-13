utilsDir="/analysisSoftware/utils_sstiefel_2024"
thisDir=`pwd`

for dir in "$utilsDir $thisDir"; do
    echo "current dir is $dir"
    if [ -n "$dir" ]; then
        trgDir="${dir}/dir_moved_compilation_files"
        mkdir $trgDir
        mv $dir/*/*/*.d */*/*.so */*/*.pcm */*/*ACLiC* ${trgDir}/.
        mv $dir/*.d *.so *.pcm *ACLiC* ${trgDir}/.
        echo "moved to ${trgDir}. ls:"
        echo `ls ${trgDir}`
    fi
done


