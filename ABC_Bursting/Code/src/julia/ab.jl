############################################# file_description ###################################
#this file contains     1. three data_scaling options
                                #  data_scaling_none
                                #  data_scaling_qPCR
                                #  data_scaling_RNAseq
#                       2. four ABC_functions
                                #  Get_single_gene_qPCR_data
                                #  Get_single_gene_RNAseq_data
                                #  Get_single_gene_RNAseq_raw_data
                                #  kolmogorov_smirnov_distance
#################################################################################################

#######################################data_scaling##############################################
function data_scaling_none(outcome::Array{Int64}, parameter::Float64)
    output = outcome[3]
    return output
end
function data_scaling_qPCR(outcome::Array{Int64}, parameter::Float64)
    if log2.(outcome[3]) > 0
        output = ((log2.(outcome[3])) * parameter)
    else
        output = 0
    end
    return output
end
function data_scaling_RNAseq(outcome::Array{Int64}, parameter::Float64)
    output = log10.((outcome[3] * parameter).+1)
    return output
end
##################################################################################################

#######################################ABC_functions##############################################
#Extracts the qPCR data for the desired gene
function Get_single_gene_qPCR_data(gene,condition,data)
    #Extract data of interest
    cell_names = data[1,2:end]
    data_values = data[2:end,1:end]
    for i in 1:69
               for j in 2:168
                          if data_values[i,j] <= 0
                                     data_values[i,j] = 0
                          end
               end
    end
    #Put data into a labelled dataframe
    dataf = convert(DataFrames.DataFrame,data_values)
    colnames = prepend!(cell_names,["Gene"])
    DataFrames.names!(dataf, [Symbol("$i") for i in colnames], allow_duplicates=true)
    #Separation of data into growth conditions
    selection = find(x->contains(x,condition),colnames)
    selection = prepend!(selection,[1])
    cond_data = dataf[:,[Symbol("$x") for x in colnames[selection]]]
    #Selection of data plot(x=test_data, Geom.histogram())
    #for comparison
    cond_data_selection = cond_data[(cond_data[:,Symbol("Gene")].==gene),:]
    test_data = Array{Float64}(convert(Array,cond_data_selection[2:end]))
    test_data = test_data[:]
    #test_data = filter(x -> x > 3,test_data)
    return test_data
end

#Extracts the RNAseq data for the desired gene and log10 transforms it
function Get_single_gene_RNAseq_data(gene,data)
    #Extract data of interest
    cell_names = data[1,1:end-1]
    data_values = data[2:end,1:end]
    #Put data into a labelled dataframe
    dataf = convert(DataFrames.DataFrame,data_values)
    colnames = prepend!(cell_names,["Gene"])
    DataFrames.names!(dataf, [Symbol("$i") for i in colnames])
    #Selection of data fro comparison
    dataf_selection =  dataf[( dataf[:,Symbol("Gene")].==gene),:]
    test_data = Array{Float64}(convert(Array, dataf_selection[2:end]))
    test_data = test_data[:]
    test_data = log10.(test_data.+1)
    return test_data
end

#Extracts the RNAseq data for the desired gene and returns raw data
function Get_single_gene_RNAseq_raw_data(gene,data)
    #Extract data of interest
    cell_names = data[1,1:end-1]
    data_values = data[2:end,1:end]
    #Put data into a labelled dataframe
    dataf = convert(DataFrames.DataFrame,data_values)
    colnames = prepend!(cell_names,["Gene"])
    DataFrames.names!(dataf, [Symbol("$i") for i in colnames])
    #Selection of data fro comparison
    dataf_selection =  dataf[( dataf[:,Symbol("Gene")].==gene),:]
    test_data = Array{Float64}(convert(Array, dataf_selection[2:end]))
    test_data = test_data[:]
    return test_data
end

#Calculates the KS distance between two distributions (approximate)
function kolmogorov_smirnov_distance(data1::Array{Float64},data2::Array{Float64})
            #Produce function which returns ecdf
            ecdf_func_1 = StatsBase.ecdf(data1)
            ecdf_func_2 = StatsBase.ecdf(data2)
            #find maximum value of both data sets for ecdf intervals
            max = maximum([data1;data2])
            intervals = max/999
            #calculate ecdf value at each interval
            ecdf_vals_1 = Array{Float64,1}(1000)
            for i in 1:1000
                        ecdf_vals_1[i]=ecdf_func_1(intervals*(i-1))
            end
            ecdf_vals_2 = Array{Float64,1}(1000)
            for i in 1:1000
                        ecdf_vals_2[i]=ecdf_func_2(intervals*(i-1))
            end
            dist = maximum(abs.(ecdf_vals_1-ecdf_vals_2))
            return dist
end
########################################################################################
