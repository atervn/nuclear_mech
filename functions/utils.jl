function zip(zipFile, srcDir)
    zdir = ZipFile.Writer(zipFile)
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

function unzip(zipFile)



end