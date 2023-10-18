using NativeFileDialog

function format_parameters()
        # Open a file dialog to select the input file
        input_file = NativeFileDialog.pick_file()
        if input_file == ""
            println("No input file selected.")
            return
        end

        # Read the contents of the input file
        contents = read(input_file, String)

        # Replace line breaks with a comma and a space
        modified_contents = replace(contents, "\r\n" => ", ")

        # Write the modified contents to the output file
        write("formatted_parameters.txt", modified_contents)

end