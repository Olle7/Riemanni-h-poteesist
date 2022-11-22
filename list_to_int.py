
def to_l(element_values,num_of_elements):
    for n in range(element_values**num_of_elements):
        for i in range(num_of_elements):
            print(  (n//element_values**(num_of_elements-1-i))%element_values   ,end=",")
        print("    ;",n)
to_l(10,4)