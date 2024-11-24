function prompt_for_integer()
    try
        print("Please enter an integer: ")
        input_str = readline()
        x = parse(Int, input_str)
        return x
    catch
        error("Invalid input. Please enter a valid integer.")
    end

end