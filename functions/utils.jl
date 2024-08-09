function zip(srcDir)
    zdir = ZipFile.Writer(srcDir*".zip")
    for (root, dirs, files) in walkdir(srcDir)
        for file in files
            filepath = joinpath(root, file)
            f = open(filepath, "r")
            content = read(f, String)
            close(f)
            zf = ZipFile.addfile(zdir, basename(filepath));
            write(zf, content)
        end
    end
    close(zdir)
end

function unzip(file)

    fileFullPath = isabspath(file) ?  file : joinpath(pwd(),file)
    basePath = dirname(fileFullPath)
    outPath = joinpath(basePath,basename(file)[1:end-4])
    isdir(outPath) ? "" : mkdir(outPath)
    zarchive = ZipFile.Reader(fileFullPath)
    for f in zarchive.files
        fullFilePath = joinpath(outPath,f.name)
        if (endswith(f.name,"/") || endswith(f.name,"\\"))
            mkdir(fullFilePath)
        else
            write(fullFilePath, read(f))
        end
    end
    close(zarchive)
end