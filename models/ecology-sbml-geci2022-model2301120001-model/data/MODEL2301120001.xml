using SymPy
# using SymEngine
using LinearAlgebra
# using CSV
# using DataFrames
using DelimitedFiles
using Plots
# using StatsPlots
# using PlotThemes
using ColorBrewer
using JLD2
# using Distributed
# using DistributedArrays
# ______________________________________________________________________________


#desired precision for numerical calculations in bits
given_precision = 120

setprecision(given_precision)

#definition of standard colours used for plotting
colour_allfemales = RGB(0/255,0/255,0/255)
colour_correlation = RGB(127/255,127/255,127/255)
colour_YLE = RGB(32/255,0/255,212/255)
colour_XShredder =  RGB(0/255,137/255,0/255)


colour_Cas9 = RGB(32/255,0/255,212/255)
colour_gRNA1 = RGB(32/255,0/255,212/255)
colour_XShredder = RGB(0/255,137/255,0/255)
colour_gRNA2 = RGB(0/255,137/255,0/255)

colour_dCas9 = RGB(32/255,0/255,212/255)
colour_dgRNA1 = RGB(32/255,0/255,212/255)
colour_dShredder = RGB(0/255,137/255,0/255)
colour_dgRNA2 = RGB(0/255,137/255,0/255)

colour_resistances = RGB(255/255,34/255,28/255)
colour_load = RGB(137/255,95/255,85/255)
colour_sexratio = RGB(255/255,0/255,255/255)



# ______________________________________________________________________________
# create a list of all genotypes in the population (in a "detailed" format)

#Y Chromosome locus 1
A_alleles = ["A","α"];
#a = wt (added later)
#A = Cas9
#α  = dysfunctional Cas9

#Y Chromosome locus 2
B_alleles = ["B","β"];
#b = wt (added later)
#B = gRNA1
#β =  dysfunctional gRNA1

#Autosome locus 1
C_alleles = ["C","ζ"];
#c = wt (added later)
#C = X-shredder
#ζ = dysfunctional X-shredder

#Autosome locus 2
D_alleles = ["D","δ"];
#d = wt (added later)
#D = gRNA2
#δ = dysfunctional gRNA2

#X Chromosome (LOCUS 1)
E_alleles = ["e", "ε", "r1"];
#e = wt editing target site
#ε = edited
#r1 = editing resistant

F_alleles = ["f", "r2"]
#f = wt shredding target site
#r2 = shredding resistant


#Define all A and B allele combinations: Creating all Y chromosome variants
Y_chromosomes = []
for A in A_alleles
    for B in B_alleles
        push!(Y_chromosomes,(A*B))
    end
end
pushfirst!(Y_chromosomes, "ab")


#Define all C and D allele combinations: Creating all autosomal chromosome variants
A_chromosomes = []
for C in C_alleles
    for D in D_alleles
        push!(A_chromosomes,(C*D))
    end
end
pushfirst!(A_chromosomes, "cd")
push!(A_chromosomes, "r3")


#Define all E and F allele combinations: Creating all X chromosome variants
X_chromosomes = []
for E in E_alleles
    for F in F_alleles
        push!(X_chromosomes,(E*F))
    end
end

#create a list of all sex genotypes
sex_chromosome_genotypes = []
for i in 1:length(Y_chromosomes)
    variant1 = Y_chromosomes[i]
    for j in 1:length(X_chromosomes)
        variant2 = X_chromosomes[j]
            push!(sex_chromosome_genotypes,([variant1 variant2]))
    end
end

for i in 1:length(X_chromosomes)
    variant1 = X_chromosomes[i]
    for j in i:length(X_chromosomes)
        variant2 = X_chromosomes[j]
        push!(sex_chromosome_genotypes, [variant1 variant2])
    end
end

#create a list of all autosomal genotypes
autosome_genotypes = []
for i in 1:length(A_chromosomes)
    variant1 = A_chromosomes[i]
    for j in i:length(A_chromosomes)
        variant2 = A_chromosomes[j]
        push!(autosome_genotypes, [variant1 variant2])
    end
end

#create a list of full detailed genotypes
genotypes_detailed = []
for i in 1:length(sex_chromosome_genotypes)
    sex_genotype = sex_chromosome_genotypes[i]
    for j in 1:length(autosome_genotypes)
        autosome_genotype = autosome_genotypes[j]
        push!(genotypes_detailed,([sex_genotype autosome_genotype]))
    end
end

# ______________________________________________________________________________
# ______________________________________________________________________________

#creating a list of all genotypes but in the simplified format

genotypes = [];

#Y Chromosome variants
Y_variants = ["Y1","Y2","Y3","Y4","Y5"];
# Y1 = fully WT Y chromosome
# Y2 = fully functional YLE (Y-linked editor) (functional Cas9 + gRNA1)
# Y3 = functional Cas9, dysfunctional gRNA1
# Y4 = dysfunctional Cas9, functional gRNA1
# Y5 = fully dysfunctional YLE

#X Chromosome variants        X
X_variants = ["X1", "X2", "X3", "X4", "X5", "X6"];
#X1 = WT
#X2 = shredding resistant
#X3 = edited
#X4 = edited, shredding resistant
#X5 = editing resistant
#X6 = editing resistant, shredding resistant

#Autosome variants
A_variants = ["A1","A2","A3","A4","A5", "A6"];
# A1 = wt
# A2 = fully functional X-shredder (functional X-shredder + gRNA2)
# A3 = functional shredder, dysfunctional gRNA2
# A4 = dysfunctional shredder, functional gRNA2
# A5 = fully dysfunctional X-shredder
# A6 = wt but homing resistant

sex_genotypes = [];

for Y in Y_variants
    for X in X_variants
        push!(sex_genotypes, [Y X])
    end
end

for i in 1:length(X_variants)
    for j in i:length(X_variants)
        push!(sex_genotypes, [X_variants[i] X_variants[j]])
    end
end

for sex_genotype in sex_genotypes
    for i in 1:length(A_variants)
        for j in i:length(A_variants)
            push!(genotypes, [sex_genotype A_variants[i] A_variants[j]])
        end
    end
end

# ______________________________________________________________________________
# ______________________________________________________________________________
# create a list of all gametes (detailed format)

sex_chromosomes = vcat(Y_chromosomes, X_chromosomes);

gametes_detailed = [];

for each_sex_chrom in sex_chromosomes
    for each_autosome in A_chromosomes
        push!(gametes_detailed, [each_sex_chrom each_autosome])
    end
end

# ______________________________________________________________________________
# ______________________________________________________________________________
# creates a list of all X chromosome variants that will be affected by shredding

shreddable_X_chromosomes = [];

for each in X_chromosomes
    if occursin("f", each)
        push!(shreddable_X_chromosomes, each)
    end
end

# ______________________________________________________________________________
# ______________________________________________________________________________

function create_selection_matrix(genotypes_detailed, )
    #calculating the fitness coefficients of the different genotypes
    #as a vector of floats (from 1 to 0) that can then be applied to the number of
    #individuals/their frequency
    #default selection coefficient is 1 for individuals that do not suffer any
    #fitness cost

    selection = sympify.(ones(Int, length(genotypes_detailed)))


    #iterating through all genotypes to check which selection costs need to
    #be applied for the specific genotypes
    for i in 1:length(genotypes_detailed)


        #first we take into account desired fitness costs
        #females homozygous for the edited target gene suffer cost s_f
        #females heterozygous for the edited target gene suffer cost s_f*h_f
        #males carrying an edited target gene on their X suffer cost s_m

        if (occursin("ε", genotypes_detailed[i][1]) && occursin("ε", genotypes_detailed[i][2]))
            selection[i] = selection[i] * (1 - Sym("s_f"))
        elseif (occursin("ε", genotypes_detailed[i][1]))
            selection[i] = selection[i] * (1 - Sym("s_f") * Sym("h_f"))
        elseif (occursin("ε", genotypes_detailed[i][2]) && (genotypes_detailed[i][1] in X_chromosomes))
            selection[i] = selection[i] * (1 - Sym("s_f") * Sym("h_f"))
        elseif (occursin("ε", genotypes_detailed[i][2]) && (genotypes_detailed[i][1] ∉ X_chromosomes))
            selection[i] = selection[i] * (1 - Sym("s_m"))
        end

        #next we take into account various forms of undesired fitness costs
        #there are expression related costs for expression of molecular components
        #and there are activity related costs that rely on activity of parts
        #they are all integrated by multiplication

        number_Cas9 = 0
        number_gRNA1 = 0
        number_shredder = 0
        number_gRNA2 = 0


        #we count how frequent the different molecular components are

        if occursin("A", genotypes_detailed[i][1])
            number_Cas9 = 1
        end

        if occursin("B", genotypes_detailed[i][1])
            number_gRNA1  = 1
        end

        if occursin("C", genotypes_detailed[i][3])
            number_shredder += 1
        end
        if occursin("C", genotypes_detailed[i][4])
            number_shredder += 1
        end

        if occursin("D", genotypes_detailed[i][3])
            number_gRNA2 += 1
        end
        if occursin("D", genotypes_detailed[i][4])
            number_gRNA2 += 1
        end

        #expression costs of proteins
        if number_Cas9 == 1
            selection[i] = selection[i] * (1 - Sym("s_a"))
        end

        if number_shredder == 2
            selection[i] = selection[i] * (1 - Sym("s_c"))*(1 - Sym("s_c"))
        elseif number_shredder == 1
            selection[i] = selection[i] * (1 - Sym("s_c"))
        end


        #expression costs of gRNAs
        if number_gRNA1 == 1
            selection[i] = selection[i] * (1 - Sym("s_b"))
        end

        if number_gRNA2 == 2
            selection[i] = selection[i] * (1 - Sym("s_b"))*(1 - Sym("s_b"))
        elseif number_gRNA2 == 1
            selection[i] = selection[i] * (1 - Sym("s_b"))
        end

        #shredder activity cost
        if number_shredder == 2
            selection[i] = selection[i] * (1 - Sym("s_e"))
        elseif number_shredder == 1
            selection[i] = selection[i] * (1 - Sym("s_e")*Sym("h_e"))
        end


        #nuclease activity costs
        if (number_Cas9 == 1 && number_gRNA1 == 1)
            selection[i] = selection[i] * (1 - Sym("s_d"))
        end

        if (number_Cas9 == 1 && number_gRNA2 >= 1)
            selection[i] = selection[i] * (1 - Sym("s_d"))
        end


    end

    return selection

end



function create_homing_matrix(genotypes_detailed, )

    homing_matrix = Matrix{Int}(I, length(genotypes_detailed), length(genotypes_detailed))
    homing_matrix = sympify.(homing_matrix)


    for i in 1:length(genotypes_detailed)

        #for homing to occur 3 conditions must be given
        #
        #(1) a functional Cas9 must be present ("A" on sex chromosome 1)
        #(2) a functional gRNA2 must be present ("D" on either autosome 1 or 2)
        #(3) one autosome needs to be wild-type
        #
        #these conditions can be met in 2 different configurations,
        #either with the wild-type autosome in position of autosome 1
        #or in autosome 2

        #BUT: genotypes are initialised to be unique and cd (wt) is at the
        #top of the genotype list, for that reason relevant cases only occur
        #with cd being at autosome 1 in the list


        Cas9 = false
        gRNA2 = false
        target = false

        if (occursin("A", genotypes_detailed[i][1]))
            Cas9 = true
        end

        if (occursin("D", genotypes_detailed[i][4]))
            gRNA2 = true
        end

        if (occursin("cd", genotypes_detailed[i][3]))
            target = true
        end



        if (Cas9 && gRNA2 && target)

            pre = deepcopy(genotypes_detailed[i])

            #determine what would be the correctly homed genotype
            post = deepcopy(pre)
            post[3] = pre[4]

            #the homing resistant genotype that occurs through mutation
            resistant = deepcopy(pre)
            resistant[3] = "r3"

            #the genotypes index numbers
            i_post = return_i(post, genotypes_detailed)
            i_resistant = return_i(resistant, genotypes_detailed)


            if (genotypes_detailed[i][4] == "CD")

                #the pre genotype will be reduced in frequency
                homing_matrix[i, i] = homing_matrix[i, i] *
                            (1 - Sym("e_h"))

                #the correctly homed genotype will increase in frequency
                homing_matrix[i, i_post] = (Sym("e_h") * (1- Sym("er_3"))*
                                        (1 - Sym("m_1"))^2)

                #resistances can occur
                homing_matrix[i, i_resistant] = Sym("e_h") * Sym("er_3")

                #accounting for the different type of mutations that might occur
                #due to copying errors
                c_mutated = deepcopy(pre)
                c_mutated[3] = "ζD"
                i_c = return_i(c_mutated, genotypes_detailed)
                homing_matrix[i, i_c] = (Sym("e_h") * (1- Sym("er_3"))*
                                        (Sym("m_1") * (1 - Sym("m_1"))))

                d_mutated = deepcopy(pre)
                d_mutated[3] = "Cδ"
                i_d = return_i(d_mutated, genotypes_detailed)
                homing_matrix[i, i_d] = (Sym("e_h") * (1- Sym("er_3"))*
                                        (Sym("m_1") * (1 - Sym("m_1"))))


                both_mutated = deepcopy(pre)
                both_mutated[3] = "ζδ"
                i_both = return_i(both_mutated, genotypes_detailed)
                homing_matrix[i, i_both] = (Sym("e_h") * (1- Sym("er_3"))*
                                                (Sym("m_1")^2))



            elseif (genotypes_detailed[i][4] == "ζD")

                homing_matrix[i, i] = homing_matrix[i, i] *
                                    (1 - Sym("e_h"))

                homing_matrix[i, i_post] = (Sym("e_h") * (1- Sym("er_3"))*
                                        (1 - Sym("m_1")))

                homing_matrix[i, i_resistant] = (Sym("e_h") *Sym("er_3"))

                d_mutated = deepcopy(pre)
                d_mutated[3] = d_mutated[4]
                d_mutated[4] = "ζδ"
                i_d = return_i(d_mutated, genotypes_detailed)
                homing_matrix[i, i_d] = (Sym("e_h") * (1- Sym("er_3"))
                                                * Sym("m_1"))

            end
        end
    end
    return homing_matrix
end


function create_recombination_matrix(genotypes_detailed)
    recombination_matrix = Matrix{Int}(I, length(genotypes_detailed), length(genotypes_detailed))
    recombination_matrix = sympify.(recombination_matrix)

    #recombination may occur between the X-linked loci (target gene + shredding target)
    #its effects only need to be considered when the genotype is changed throug
    #recombination
    #for example for homozygous individuals this would not be the case
    #only double heterozygous individuals are of relevance

    for i in 1:length(genotypes_detailed)

        if (genotypes_detailed[i][1] in X_chromosomes)

            A = []
            B = []

            if occursin("e",genotypes_detailed[i][1])
                push!(A, "e")
            elseif occursin("ε",genotypes_detailed[i][1])
                push!(A, "ε")
            elseif occursin("r1",genotypes_detailed[i][1])
                push!(A, "r1")
            end

            if occursin("e",genotypes_detailed[i][2])
                push!(B, "e")
            elseif occursin("ε",genotypes_detailed[i][2])
                push!(B, "ε")
            elseif occursin("r1",genotypes_detailed[i][2])
                push!(B, "r1")
            end

            if (A != B)

                if occursin("f",genotypes_detailed[i][1])
                    push!(A, "f")
                elseif occursin("r2",genotypes_detailed[i][1])
                    push!(A, "r2")
                end

                if occursin("f",genotypes_detailed[i][2])
                    push!(B, "f")
                elseif occursin("r2",genotypes_detailed[i][2])
                    push!(B, "r2")
                end

                if (A[2] != B[2])
                    post = deepcopy(genotypes_detailed[i])
                    post[1] = string(A[1], B[2])
                    post[2] = string(B[1], A[2])
                    i_post = return_i(post, genotypes_detailed)

                    # println(genotypes_detailed[i][1], " and ", genotypes_detailed[i][2])
                    # println("become")
                    # println(post[1], " and ", post[2])

                    recombination_matrix[i,i] = (1 - Sym("r"))
                    recombination_matrix[i,i_post] = Sym("r")

                end

            end


        end


    end

    return recombination_matrix
end


function create_editing_matrix(genotypes_detailed, )

    editing_matrix = Matrix{Int}(I, length(genotypes_detailed), length(genotypes_detailed))
    editing_matrix = sympify.(editing_matrix)


    for i in 1:length(genotypes_detailed)

        # for editing to occur three conditions must be fulfilled:
        #(1) a functional Cas9 must be present
        #(2) a functional gRNA1 must be present
        #(3) an editable X chromosome must be present

        if ((genotypes_detailed[i][1] == "AB") && (occursin("e",genotypes_detailed[i][2])))
            post = deepcopy(genotypes_detailed[i])
            post[2] = replace(post[2], "e" => "ε")

            resistant = deepcopy(genotypes_detailed[i])
            resistant[2] = replace(resistant[2], "e" => "r1")

            i_post = return_i(post, genotypes_detailed)
            i_resistant = return_i(resistant, genotypes_detailed)

            editing_matrix[i,i] = (1 - Sym("e_e"))

            editing_matrix[i, i_post] = (Sym("e_e") * (1- Sym("er_1")))

            editing_matrix[i, i_resistant] = (Sym("e_e")* Sym("er_1"))


        end
    end
    return editing_matrix
end



function create_mutation_matrix(genotypes_detailed, )

    #this is the matrix describing the effect of random background mutations
    #we assume that may occur anywhere in the genome and is the same for all elements

    mutation_matrix = Matrix{Int}(I, length(genotypes_detailed), length(genotypes_detailed))
    mutation_matrix = sympify.(mutation_matrix)

    for i in 1:length(genotypes_detailed)

        sex1 = [genotypes_detailed[i][1] 0]
        sex2 = [genotypes_detailed[i][2] 0]
        aut1 = [genotypes_detailed[i][3] 0]
        aut2 = [genotypes_detailed[i][4] 0]

        #for every genotype it is determined which functional elements exist that
        #can become dysfunctional as a consequence of mutation
        #further, it is saved how many mutations this would require to occur at once

        if (genotypes_detailed[i][1] == "AB")
            sex1 = vcat(sex1, ["αB" 1; "Aβ" 1; "αβ" 2])
        elseif  (occursin("A", genotypes_detailed[i][1]))
            mutation = replace(genotypes_detailed[i][1], "A" => "α")
            sex1 = vcat(sex1, [mutation 1])
        elseif  (occursin("B", genotypes_detailed[i][1]))
            mutation = replace(genotypes_detailed[i][1], "B" => "β")
            sex1 = vcat(sex1, [mutation 1])
        end

        if (genotypes_detailed[i][3] == "CD")
            aut1 = vcat(aut1, ["ζD" 1; "Cδ" 1; "ζδ" 2])
        elseif  (occursin("C", genotypes_detailed[i][3]))
            mutation = replace(genotypes_detailed[i][3], "C" => "ζ")
            aut1 = vcat(aut1, [mutation 1])
        elseif  (occursin("D", genotypes_detailed[i][3]))
            mutation = replace(genotypes_detailed[i][3], "D" => "δ")
            aut1 = vcat(aut1, [mutation 1])
        end

        if (genotypes_detailed[i][4] == "CD")
            aut2 = vcat(aut2, ["ζD" 1; "Cδ" 1; "ζδ" 2])
        elseif  (occursin("C", genotypes_detailed[i][4]))
            mutation = replace(genotypes_detailed[i][4], "C" => "ζ")
            aut2 = vcat(aut2, [mutation 1])
        elseif  (occursin("D", genotypes_detailed[i][4]))
            mutation = replace(genotypes_detailed[i][4], "D" => "δ")
            aut2 = vcat(aut2, [mutation 1])
        end


        mutated_genotypes = Any[0 0]
        mutated_genotypes[1] = genotypes_detailed[i]

        #finally all combinations of mutations in the individual elements are
        #taken to together to determine all possible genotypes that may arise
        #throug mutation and also how many  mutations would have to
        #occur simultaneously to give the corresponding genotype in a single step

        for j in 1:length(sex1[:,1])
            for k in 1:length(sex2[:,1])
                for l in 1:length(aut1[:,1])
                    for m in 1:length(aut2[:,1])
                        number_of_mutations = sex1[j,2] + sex2[k,2] + aut1[l,2] + aut2[m,2]
                        mutation = [sex1[j,1] sex2[k,1] aut1[l,1] aut2[m,1]]

                        mutated_genotype = Any[0 number_of_mutations]
                        mutated_genotype[1] = mutation

                        if (number_of_mutations > 0)
                            mutated_genotypes = vcat(mutated_genotypes, mutated_genotype)
                        end
                    end
                end
            end
        end

        #next, go through all genotypes and for each modify mutation matrix
        #values for post and pre mutation genotypes accordingly

        no_mutatable_elements = mutated_genotypes[end,2]

        for n in 1:length(mutated_genotypes[:,1])

            if (mutated_genotypes[n,2] > 0)

                i_post = return_i(mutated_genotypes[n,1], genotypes_detailed)
                mutation_matrix[i,i] = (1 - Sym("m_2"))^(mutated_genotypes[n,2])

                no_not_mutated = no_mutatable_elements - mutated_genotypes[n,2]

                if no_not_mutated > 0
                    mutation_matrix[i, i_post] += (Sym("m_2"))^(mutated_genotypes[n,2]) *
                                                 (1 - Sym("m_2"))^(no_not_mutated)
                else
                    mutation_matrix[i, i_post] += (Sym("m_2"))^(mutated_genotypes[n,2])
                end
            end

        end

    end

    #simplifying terms if possible
    for a in length(mutation_matrix[:,1])
        for b in length(mutation_matrix[1,:])
            mutation_matrix[a,b] = simplify(mutation_matrix[a,b])
        end
    end


    return mutation_matrix
end


function create_gamete_matrix(genotypes_detailed, gametes_detailed, )
    gamete_matrix = zeros(Float64,length(genotypes_detailed), length(gametes_detailed))
    gamete_matrix = sympify.(gamete_matrix)


    for i in 1:length(genotypes_detailed)
        genotype = genotypes_detailed[i]
        gametes = [[[genotype[1] genotype[3]]]
                    [[genotype[1] genotype[4]]]
                    [[genotype[2] genotype[3]]]
                    [[genotype[2] genotype[4]]]
                        ]

        if (occursin("C", genotype[3]) && occursin("C", genotype[4]))
            shredder = 2
        elseif (occursin("C", genotype[3]) || occursin("C", genotype[4]))
            shredder = 1
        else
            shredder = 0
        end

        if (genotype[1] in X_chromosomes)
            female = true
        else
            female = false
        end

        if (genotype[2] in shreddable_X_chromosomes)
            shreddable = true
        else
            shreddable = false
        end


        if (female || (shredder == 0) || ! shreddable )

            for gamete in gametes
                i_gamete = return_i(gamete, gametes_detailed)
                if (gamete_matrix[i, i_gamete] == 0)
                    frequency = count(x -> (x == gamete), gametes)
                    gamete_matrix[i, i_gamete] = frequency / 4
                end
            end

        else
            resistant1 = deepcopy(gametes[3])
            resistant2 = deepcopy(gametes[4])

            resistant1[1] = replace(resistant1[1], "f" => "r2")
            resistant2[1] = replace(resistant2[1], "f" => "r2")

            i_resistant1 = return_i(resistant1, gametes_detailed)
            i_resistant2 = return_i(resistant2, gametes_detailed)

            i_1 = return_i(gametes[1], gametes_detailed)
            i_2 = return_i(gametes[2], gametes_detailed)
            i_3 = return_i(gametes[3], gametes_detailed)
            i_4 = return_i(gametes[4], gametes_detailed)


            if (shredder == 2) && (genotype[1] == "AB")

                gamete_matrix[i, i_1] += (1/(4 - 2*Sym("e_s")*(1-Sym("er_2"))))
                gamete_matrix[i, i_2] += (1/(4 - 2*Sym("e_s")*(1-Sym("er_2"))))
                gamete_matrix[i, i_3] += (1 - Sym("e_s"))/(4 - 2*Sym("e_s")*(1-Sym("er_2")))
                gamete_matrix[i, i_4] +=  (1 - Sym("e_s"))/(4 - 2*Sym("e_s")*(1-Sym("er_2")))
                gamete_matrix[i, i_resistant1] +=  (Sym("e_s")*Sym("er_2"))/(4 - 2*Sym("e_s")*(1-Sym("er_2")))
                gamete_matrix[i, i_resistant2] +=  (Sym("e_s")*Sym("er_2"))/(4 - 2*Sym("e_s")*(1-Sym("er_2")))


            elseif (shredder == 2)  && (genotype[1] !== "AB")

                gamete_matrix[i, i_1] += (1/(4 - 2*Sym("e_s")*Sym("c")*(1-Sym("er_2"))))
                gamete_matrix[i, i_2] += (1/(4 - 2*Sym("e_s")*Sym("c")*(1-Sym("er_2"))))
                gamete_matrix[i, i_3] += (1 - Sym("e_s")*Sym("c"))/(4 - 2*Sym("e_s")*Sym("c")*(1-Sym("er_2")))
                gamete_matrix[i, i_4] +=  (1 - Sym("e_s")*Sym("c"))/(4 - 2*Sym("e_s")*Sym("c")*(1-Sym("er_2")))
                gamete_matrix[i, i_resistant1] +=  (Sym("e_s")*Sym("er_2")*Sym("c"))/(4 - 2*Sym("e_s")*Sym("c")*(1-Sym("er_2")))
                gamete_matrix[i, i_resistant2] +=  (Sym("e_s")*Sym("er_2")*Sym("c"))/(4 - 2*Sym("e_s")*Sym("c")*(1-Sym("er_2")))


            elseif (shredder == 1 ) && (genotype[1] == "AB")

                gamete_matrix[i, i_1] += (1/(4 - 2*Sym("e_s")*Sym("h_e2")*(1-Sym("er_2"))))
                gamete_matrix[i, i_2] += (1/(4 - 2*Sym("e_s")*Sym("h_e2")*(1-Sym("er_2"))))
                gamete_matrix[i, i_3] += (1 - Sym("e_s")*Sym("h_e2"))/(4 - 2*Sym("e_s")*Sym("h_e2")*(1-Sym("er_2")))
                gamete_matrix[i, i_4] +=  (1 - Sym("e_s")*Sym("h_e2"))/(4 - 2*Sym("e_s")*Sym("h_e2")*(1-Sym("er_2")))
                gamete_matrix[i, i_resistant1] +=  Sym("e_s")*Sym("h_e2")*Sym("er_2")/(4 - 2*Sym("e_s")*Sym("h_e2")*(1-Sym("er_2")))
                gamete_matrix[i, i_resistant2] +=  Sym("e_s")*Sym("h_e2")*Sym("er_2")/(4 - 2*Sym("e_s")*Sym("h_e2")*(1-Sym("er_2")))



            elseif (shredder == 1)  && (genotype[1] !== "AB")

                gamete_matrix[i, i_1] += (1/(4 - 2*Sym("e_s")*Sym("h_e2")*Sym("c")*(1-Sym("er_2"))))
                gamete_matrix[i, i_2] += (1/(4 - 2*Sym("e_s")*Sym("h_e2")*Sym("c")*(1-Sym("er_2"))))
                gamete_matrix[i, i_3] += (1 - Sym("e_s")*Sym("h_e2")*Sym("c"))/(4 - 2*Sym("e_s")*Sym("h_e2")*Sym("c")*(1-Sym("er_2")))
                gamete_matrix[i, i_4] +=  (1 - Sym("e_s")*Sym("h_e2")*Sym("c"))/(4 - 2*Sym("e_s")*Sym("h_e2")*Sym("c")*(1-Sym("er_2")))

                # println(i_1)
                # println(i_2)
                # println(i_3)
                # println(i_4)
                # println(i_resistant1)
                # println(i_resistant2)

                gamete_matrix[i, i_resistant1] +=  (Sym("e_s")*Sym("h_e2")*Sym("c")*Sym("er_2"))/(4 - 2*Sym("e_s")*Sym("h_e2")*Sym("c")*(1-Sym("er_2")))
                gamete_matrix[i, i_resistant2] +=  (Sym("e_s")*Sym("h_e2")*Sym("c")*Sym("er_2"))/(4 - 2*Sym("e_s")*Sym("h_e2")*Sym("c")*(1-Sym("er_2")))



            end
        end
    end

    #     elseif ((! female) && (shredder == 2))
    #
    #         for gamete in gametes
    #             i_gamete = return_i(gamete, gametes_detailed)
    #             frequency = count(x -> (x == gamete), gametes)
    #
    #             if ((gamete_matrix[i, i_gamete] == 0)
    #                 && (gamete[1] in X_chromosomes) && (genotypes_detailed[i][1] == "AB"))
    #                 gamete_matrix[i, i_gamete] = (frequency / 4) * (1 - Sym("e_s"))
    #
    #             elseif ((gamete_matrix[i, i_gamete]) == 0
    #                 && (gamete[1] ∉ X_chromosomes) && (genotypes_detailed[i][1] == "AB"))
    #                 gamete_matrix[i, i_gamete] = (frequency / 4) * (1 + Sym("e_s"))
    #
    #             elseif ((gamete_matrix[i, i_gamete] == 0)
    #                     && (gamete[1] in X_chromosomes) && (genotypes_detailed[i][1] != "AB"))
    #                     gamete_matrix[i, i_gamete] = (frequency / 4) *
    #                             (1 - Sym("e_s")*Sym("c"))
    #
    #             elseif ((gamete_matrix[i, i_gamete]) == 0
    #                 && (gamete[1] ∉ X_chromosomes) && (genotypes_detailed[i][1] != "AB"))
    #
    #                 gamete_matrix[i, i_gamete] = (frequency / 4) *
    #                             (1 + Sym("e_s") * Sym("c"))
    #
    #             end
    #         end
    #
    #     elseif ((! female) && (shredder == 1))
    #
    #         for gamete in gametes
    #             i_gamete = return_i(gamete, gametes_detailed)
    #             frequency = count(x -> (x == gamete), gametes)
    #             # println(gamete, " in ", gametes)
    #             # println(frequency)
    #             if (gamete_matrix[i, i_gamete] == 0 &&
    #                  (gamete[1] in X_chromosomes) && (genotypes_detailed[i][1] == "AB"))
    #
    #                 gamete_matrix[i, i_gamete] = (frequency / 4) *
    #                             (1 - Sym("e_s")*Sym("h_e2"))
    #
    #             elseif (gamete_matrix[i, i_gamete] == 0 &&
    #                 (gamete[1] ∉ X_chromosomes) && (genotypes_detailed[i][1] == "AB"))
    #
    #                 gamete_matrix[i, i_gamete] = (frequency / 4) *
    #                         (1 + Sym("e_s")*Sym("h_e2"))
    #
    #             elseif (gamete_matrix[i, i_gamete] == 0 &&
    #                         (gamete[1] in X_chromosomes) && (genotypes_detailed[i][1] != "AB"))
    #
    #                     gamete_matrix[i, i_gamete] = (frequency / 4) *
    #                             (1 - Sym("e_s")*Sym("h_e2")*Sym("c"))
    #
    #             elseif (gamete_matrix[i, i_gamete] == 0 &&
    #                     (gamete[1] ∉ X_chromosomes) && (genotypes_detailed[i][1] != "AB"))
    #
    #                     gamete_matrix[i, i_gamete] = (frequency / 4) *
    #                         (1 + Sym("e_s")*Sym("h_e2")*Sym("c"))
    #             end
    #         end
    #     end
    # end


    break_point = 0

    for i in 1:length(genotypes_detailed)

        if ((break_point == 0) && (genotypes_detailed[i][1] in X_chromosomes))
            break_point = i
            break
        end
    end


    empty_matrix_up = zeros(Int,(break_point - 1), length(gametes_detailed))

    empty_matrix_down = zeros(Int,(length(genotypes_detailed[break_point:end])), length(gametes_detailed))

    sperm_matrix = vcat(gamete_matrix[1:(break_point - 1),:], empty_matrix_down)
    egg_matrix = vcat(empty_matrix_up ,gamete_matrix[break_point:end,:])

    return sperm_matrix, egg_matrix

end


function create_zygote_vector(sperm_vector, egg_vector)

    zygote_matrix = zeros(BigFloat, length(genotypes_detailed))

    for i in 1:length(sperm_vector)
        if (sperm_vector[i] != 0.0)
            for k in 1:length(egg_vector)
                if (egg_vector[k] != 0.0)
                    genotype = [gametes_detailed[i][1] gametes_detailed[k][1] gametes_detailed[i][2] gametes_detailed[k][2]]
                    #println(genotype)
                    freq = sperm_vector[i] * egg_vector[k]
                    i_genotype = return_i(genotype, genotypes_detailed)
                    #println(i_genotype)
                    zygote_matrix[i_genotype] += freq
                end
            end
        end
    end

    return zygote_matrix

end


function create_zygote_vector_sym(sperm_vector, egg_vector)

    zygote_matrix = zeros(Float64, length(genotypes_detailed))
    zygote_matrix = sympify.(zygote_matrix)

    for i in 1:length(sperm_vector)
        if (sperm_vector[i] != 0.0)
            for k in 1:length(egg_vector)
                if (egg_vector[k] != 0.0)
                    genotype = [gametes_detailed[i][1] gametes_detailed[k][1] gametes_detailed[i][2] gametes_detailed[k][2]]
                    #println(genotype)
                    freq = sperm_vector[i] * egg_vector[k]
                    i_genotype = return_i(genotype, genotypes_detailed)
                    #println(i_genotype)
                    zygote_matrix[i_genotype] += freq
                end
            end
        end
    end

    return zygote_matrix

end




function compare_genotypes(geno1, geno2)
    #this function takes two genotypes (represented as the arrays of strings)
    #it compares them and returns TRUE if they are the same

    result = true
    for each in geno1
        if each ∉ geno2
            result = false
            break
        end
    end

    for each in geno2
        if each ∉ geno1
            result = false
            break
        end
    end
    return result
end


function return_i(geno1, genotypes)

    #this function takes a genotype and a list of genotypes
    #and it returns the (first) index number of the genotype in the list
    #this allows it to easily find a genotype's index number
    #implicitly assumes that the list only contains unique genotypes (which is the case)

    for i in 1:length(genotypes)
        if (compare_genotypes(geno1, genotypes[i]))
            return i
        end
    end
end

function multiplication(a,b)
    #this function performs the multiplication of genotype vector with
    #the tables that describe the different processes
    #it takes the genotype vector "b" and the table/matrix "a"

    c = zeros(length(a[1,:]))

    for i in 1:length(a[:,1])
        for k in 1:length(a[1,:])
            if (a[i,k] != 0) && (b[i] != 0)
                c[k] = c[k] + a[i,k]*b[i]
            end
        end
    end

    return c
end

function multiplication_sym(a,b)
    #this is the same as multiplication() but it returns results as type of Symbol

    c = zeros(length(a[1,:]))
    c = sympify.(c)

    for i in 1:length(a[:,1])
        for k in 1:length(a[1,:])
            if (a[i,k] != 0) && (b[i] != 0)
                c[k] = c[k] + a[i,k]*b[i]
            end
        end
    end

    return c
end


function multiplication_optimised(a,b,helper)
    # this is a multiplication function that is more efficient than multiplication()
    # because it uses a before-hand calculated "helper table" that defines
    # which entries in each table need to be multiplied, i.e. are not zero

    output_vector = BigFloat.(zeros(length(a[1,:])))

    for i in 1:length(a[:,1])
        for k in helper[i]
                output_vector[k] = output_vector[k] + a[i,k]*b[i]
        end
    end

    return output_vector
end


function plausibility_check(args ...)
    # this function performs a simple plausibility control on tables by checking
    # whether all rows add up to 1 or 0
    # i.e. whether the tables are consistent and dont lead to gain or loss of
    # individuals

    for arg in args
        ok = true
        for i in 1:length(arg[:,1])
            a = sum(arg[i,:])
            if typeof(a) == Sym
                a = simplify(a)
            end
            #a = round(sum(arg[i,:]), digits = 15)
            if (a ∉ [1.0  0.0 1 0])
                ok = false
                println(a, " in index ", i)
                break
            end
        end
    println("Plausibility check ok? ", ok)
    end
end



function plausibility_check2(args ...)
    # this function is another plausibility check function but more detailed
    # and it returns information in which rows a wrong value occured

    for arg in args
        ok = true
        for i in 1:length(arg[:,1])
            a = sum(arg[i,:])
            if (a ∉ [1.0  0.0])
                ok = false
                println(a, " in index ", i)
            end

            for j in 1:length(arg[1,:])
                if (arg[i,j] > 1) || (arg[i,j] < 0)
                    ok = false
                    println(arg[i,j], " is not between 0 and 1")
                end
            end
        end
    println("Plausibility check ok? ", ok)
    end
end


# this is the generic function for calculating time courses of a
# provided population and time
# as a single, independent population
# it also calculates in every generation what the reproductive output would be
# if there was no editing, homing, fitness effects etc to later calculate genetic load
function timecourse(t, genotype_vector; initial_size = 1.0, number_of_releases = 1,
    release_genotype = [], release_size = 0 )
    global current_Parameters, current_matrices, wildtype_matrices

    store_genotypes = deepcopy(genotype_vector)
    store_zygotes = BigFloat.(zeros(length(genotype_vector)))
    store_eggs = BigFloat.(zeros(length(gametes_detailed)))
    store_sperm = BigFloat.(zeros(length(gametes_detailed)))
    store_wildtype = deepcopy(genotype_vector)

    Rm = current_Parameters["Rm"]
    theta = current_Parameters["theta"]

    #number of eggs ("f")
    f = (Rm * 2.0) / theta

    #alpha is a constant determining the density dependent mortality (half maximal)
    alpha = initial_size * f / (Rm - 1.0)

    for i in 1:t
        genotype_vector_wildtype = deepcopy(genotype_vector)

        if (number_of_releases > 1 && i <= number_of_releases)
            genotype_vector[return_i(release_genotype, genotypes_detailed)] += release_size
        end

        genotype_vector = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector, mutation_matrix_helper)
        genotype_vector_wildtype = multiplication_optimised(wildtype_matrices["mutation_matrix"], genotype_vector_wildtype, mutation_matrix_helper)

        genotype_vector = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector, homing_matrix_helper)
        genotype_vector_wildtype = multiplication_optimised(wildtype_matrices["homing_matrix"], genotype_vector_wildtype, homing_matrix_helper)

        genotype_vector = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector, editing_matrix_helper)
        genotype_vector_wildtype = multiplication_optimised(wildtype_matrices["editing_matrix"], genotype_vector_wildtype, editing_matrix_helper)

        genotype_vector = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector, recombination_matrix_helper)
        genotype_vector_wildtype = multiplication_optimised(wildtype_matrices["recombination_matrix"], genotype_vector_wildtype, recombination_matrix_helper)

        sperm = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector, sperm_matrix_helper)
        sperm_wildtype = multiplication_optimised(wildtype_matrices["sperm_matrix"], genotype_vector_wildtype, sperm_matrix_helper)

        total_sperm = sum(sperm)
        total_sperm_wildtype = sum(sperm_wildtype)

        sperm = sperm./total_sperm
        sperm_wildtype = sperm_wildtype./total_sperm_wildtype

        eggs = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector, egg_matrix_helper)
        eggs_wildtype =multiplication_optimised(wildtype_matrices["egg_matrix"], genotype_vector_wildtype, egg_matrix_helper)

        eggs = eggs .* f
        eggs_wildtype = eggs_wildtype .* f

        store_sperm = hcat(store_sperm, sperm)
        store_eggs = hcat(store_eggs, eggs)

        genotype_vector = create_zygote_vector(sperm,eggs)
        genotype_vector_wildtype = create_zygote_vector(sperm_wildtype, eggs_wildtype)

        store_zygotes = hcat(store_zygotes, genotype_vector)

        #println("RUN ", i,)
        #println(sperm)
        #println(eggs)

        total_zygotes = sum(genotype_vector)
        total_zygotes_wildtype = sum(genotype_vector_wildtype)

        genotype_vector = genotype_vector .* (theta * alpha / (alpha + total_zygotes))
        genotype_vector_wildtype = genotype_vector_wildtype .* (theta * alpha / (alpha + total_zygotes_wildtype))

        genotype_vector = current_matrices["selection_matrix"] .* genotype_vector
        genotype_vector_wildtype = wildtype_matrices["selection_matrix"] .* genotype_vector_wildtype

        store_genotypes = hcat(store_genotypes, genotype_vector)
        store_wildtype = hcat(store_wildtype, genotype_vector_wildtype)

    end

    run = Dict(
        "genotypes" => store_genotypes,
        "zygotes" => store_zygotes,
        "eggs" => store_eggs,
        "sperm" => store_sperm,
        "wildtype" => store_wildtype,
    )

    return run
end



function timecourse_parameterscan(t, genotype_vector; initial_size = 1.0, threshold = 0.01)
    global current_Parameters, current_matrices

    max_suppression = 1
    duration_of_protection = 0
    start_duration = 0

    Rm = current_Parameters["Rm"]
    theta = current_Parameters["theta"]

    #number of eggs ("f")
    f = (Rm * 2.0) / theta

    #alpha is a constant determining the density dependent mortality (half maximal)
    alpha = initial_size * f / (Rm - 1.0)

    sex_breakpoint = 0
    for m in 1:length(genotypes_detailed)
        if genotypes_detailed[m][1] in X_chromosomes
            sex_breakpoint = m
            break
        end
    end

    # println(sex_breakpoint)

    for i in 1:t

        genotype_vector = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector, mutation_matrix_helper)
        genotype_vector = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector, homing_matrix_helper)
        genotype_vector = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector, editing_matrix_helper)
        genotype_vector = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector, recombination_matrix_helper)

        sperm = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector, sperm_matrix_helper)
        total_sperm = sum(sperm)
        sperm = sperm./total_sperm

        eggs = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector, egg_matrix_helper)
        eggs = eggs .* f

        # println(sum(sperm))
        # println(sum(eggs))

        genotype_vector = create_zygote_vector(sperm,eggs)

        total_zygotes = sum(genotype_vector)

        genotype_vector = genotype_vector .* (theta * alpha / (alpha + total_zygotes))

        genotype_vector = current_matrices["selection_matrix"] .* genotype_vector

        # println(sum(genotype_vector))
        allfemales = 0
        for j in sex_breakpoint:length(genotypes_detailed)
            allfemales = allfemales + genotype_vector[j]
        end


        if allfemales < max_suppression
            max_suppression = allfemales
        end

        # println(allfemales)

        if (start_duration == 0) && (allfemales < threshold)
            start_duration = i
        elseif (start_duration != 0) && (allfemales > threshold)
            duration_of_protection = i - start_duration
            break
        end
    end

    if (start_duration != 0) && (duration_of_protection == 0)
        duration_of_protection = t - start_duration
    end


    return max_suppression, duration_of_protection
end


function timecourse_parameterscan_comprehensive(t, genotype_vector; initial_size = 1.0, threshold = 0.01)
    global current_Parameters, current_matrices

    max_suppression = BigFloat("1")
    duration_of_protection = 0
    start_duration = 0
    transgenic_Y = BigFloat("0")
    transgenic_S = BigFloat("0")

    Rm = current_Parameters["Rm"]
    theta = current_Parameters["theta"]

    #number of eggs ("f")
    f = (Rm * 2.0) / theta

    #alpha is a constant determining the density dependent mortality (half maximal)
    alpha = initial_size * f / (Rm - 1.0)

    # sex_breakpoint = 0
    # for m in 1:length(genotypes_detailed)
    #     if genotypes_detailed[m][1] in X_chromosomes
    #         sex_breakpoint = m
    #         break
    #     end
    # end

    # println(sex_breakpoint)

    for i in 1:t

        genotype_vector = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector, mutation_matrix_helper)
        genotype_vector = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector, homing_matrix_helper)
        genotype_vector = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector, editing_matrix_helper)
        genotype_vector = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector, recombination_matrix_helper)

        sperm = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector, sperm_matrix_helper)
        total_sperm = sum(sperm)
        sperm = sperm./total_sperm

        eggs = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector, egg_matrix_helper)
        eggs = eggs .* f

        # println(sum(sperm))
        # println(sum(eggs))

        genotype_vector = create_zygote_vector(sperm,eggs)

        total_zygotes = sum(genotype_vector)

        genotype_vector = genotype_vector .* (theta * alpha / (alpha + total_zygotes))

        genotype_vector = current_matrices["selection_matrix"] .* genotype_vector

        # println(sum(genotype_vector))
        allfemales = BigFloat("0")
        allmales = BigFloat("0")
        cur_transgenic_Y = BigFloat("0")
        cur_transgenic_S = BigFloat("0")


        for k in 1:length(genotype_vector)

            if genotypes_detailed[k][1] in X_chromosomes
                allfemales += genotype_vector[k]
            else
                allmales += genotype_vector[k]
            end

            if occursin("AB", genotypes_detailed[k][1])
                cur_transgenic_Y +=  genotype_vector[k]
            elseif occursin("A", genotypes_detailed[k][1])
                cur_transgenic_Y +=  genotype_vector[k]
            elseif occursin("B", genotypes_detailed[k][1])
                cur_transgenic_Y +=  genotype_vector[k]
            end


            if ((occursin("C", genotypes_detailed[k][3]) || occursin("D", genotypes_detailed[k][3]))
                &&  (occursin("C", genotypes_detailed[k][4])) || occursin("D", genotypes_detailed[k][4]))
                cur_transgenic_S += 2 * genotype_vector[k]

            elseif ((occursin("C", genotypes_detailed[k][3]) || occursin("D", genotypes_detailed[k][3]))
                ||  (occursin("C", genotypes_detailed[k][4])) || occursin("D", genotypes_detailed[k][4]))
                cur_transgenic_S +=  genotype_vector[k]
            end

        end

        cur_transgenic_Y = cur_transgenic_Y / allmales
        cur_transgenic_S = cur_transgenic_S / (allmales * 2 + allfemales * 2)


        if allfemales < max_suppression
            max_suppression = allfemales
        end

        if cur_transgenic_Y > transgenic_Y
            transgenic_Y = cur_transgenic_Y
        end

        if cur_transgenic_S > transgenic_S
            transgenic_S = cur_transgenic_S
        end


        # println(allfemales)

        if (start_duration == 0) && (allfemales < threshold)
            start_duration = i
        elseif (start_duration != 0) && (allfemales > threshold)
            duration_of_protection = i - start_duration
            break
        end
    end

    if (start_duration != 0) && (duration_of_protection == 0)
        duration_of_protection = t - start_duration
    end


    return max_suppression, duration_of_protection, transgenic_Y, transgenic_S
end


function timecourse_parameterscan_preexisting_resistance(t, genotype_vector; initial_size = 1.0, threshold = 0.01)
    global current_Parameters, current_matrices

    max_suppression = BigFloat("1")
    t_max_suppression = 0
    t_1e_2 = 0
    t_1e_3 = 0
    t_1e_4 = 0
    t_1e_5 = 0
    t_1e_6 = 0
    t_1e_7 = 0
    t_1e_8 = 0
    t_1e_9 = 0
    t_1e_10 = 0
    duration_of_protection = 0
    start_duration = 0
    transgenic_Y = BigFloat("0")
    transgenic_S = BigFloat("0")

    Rm = current_Parameters["Rm"]
    theta = current_Parameters["theta"]

    #number of eggs ("f")
    f = (Rm * 2.0) / theta

    #alpha is a constant determining the density dependent mortality (half maximal)
    alpha = initial_size * f / (Rm - 1.0)

    # sex_breakpoint = 0
    # for m in 1:length(genotypes_detailed)
    #     if genotypes_detailed[m][1] in X_chromosomes
    #         sex_breakpoint = m
    #         break
    #     end
    # end

    # println(sex_breakpoint)

    for i in 1:t

        genotype_vector = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector, mutation_matrix_helper)
        genotype_vector = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector, homing_matrix_helper)
        genotype_vector = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector, editing_matrix_helper)
        genotype_vector = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector, recombination_matrix_helper)

        sperm = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector, sperm_matrix_helper)
        total_sperm = sum(sperm)
        sperm = sperm./total_sperm

        eggs = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector, egg_matrix_helper)
        eggs = eggs .* f

        # println(sum(sperm))
        # println(sum(eggs))

        genotype_vector = create_zygote_vector(sperm,eggs)

        total_zygotes = sum(genotype_vector)

        genotype_vector = genotype_vector .* (theta * alpha / (alpha + total_zygotes))

        genotype_vector = current_matrices["selection_matrix"] .* genotype_vector

        # println(sum(genotype_vector))
        allfemales = BigFloat("0")
        allmales = BigFloat("0")
        cur_transgenic_Y = BigFloat("0")
        cur_transgenic_S = BigFloat("0")


        for k in 1:length(genotype_vector)

            if genotypes_detailed[k][1] in X_chromosomes
                allfemales += genotype_vector[k]
            else
                allmales += genotype_vector[k]
            end

            if occursin("AB", genotypes_detailed[k][1])
                cur_transgenic_Y +=  genotype_vector[k]
            elseif occursin("A", genotypes_detailed[k][1])
                cur_transgenic_Y +=  genotype_vector[k]
            elseif occursin("B", genotypes_detailed[k][1])
                cur_transgenic_Y +=  genotype_vector[k]
            end


            if ((occursin("C", genotypes_detailed[k][3]) || occursin("D", genotypes_detailed[k][3]))
                &&  (occursin("C", genotypes_detailed[k][4])) || occursin("D", genotypes_detailed[k][4]))
                cur_transgenic_S += 2 * genotype_vector[k]

            elseif ((occursin("C", genotypes_detailed[k][3]) || occursin("D", genotypes_detailed[k][3]))
                ||  (occursin("C", genotypes_detailed[k][4])) || occursin("D", genotypes_detailed[k][4]))
                cur_transgenic_S +=  genotype_vector[k]
            end

        end

        cur_transgenic_Y = cur_transgenic_Y / allmales
        cur_transgenic_S = cur_transgenic_S / (allmales * 2 + allfemales * 2)


        if allfemales < max_suppression
            max_suppression = allfemales
            t_max_suppression = i
        end

        if (allfemales < 1e-2) && (t_1e_2 == 0)
            t_1e_2 = i
        end

        if (allfemales < 1e-3) && (t_1e_3 == 0)
            t_1e_3 = i
        end

        if (allfemales < 1e-4) && (t_1e_4 == 0)
            t_1e_4 = i
        end

        if (allfemales < 1e-5) && (t_1e_5 == 0)
            t_1e_5 = i
        end

        if (allfemales < 1e-6) && (t_1e_6 == 0)
            t_1e_6 = i
        end

        if (allfemales < 1e-7) && (t_1e_7 == 0)
            t_1e_7 = i
        end

        if (allfemales < 1e-8) && (t_1e_8 == 0)
            t_1e_8 = i
        end

        if (allfemales < 1e-9) && (t_1e_9 == 0)
            t_1e_9 = i
        end

        if (allfemales < 1e-10) && (t_1e_10 == 0)
            t_1e_10 = i
        end

        if cur_transgenic_Y > transgenic_Y
            transgenic_Y = cur_transgenic_Y
        end

        if cur_transgenic_S > transgenic_S
            transgenic_S = cur_transgenic_S
        end


        # println(allfemales)

        if (start_duration == 0) && (allfemales < threshold)
            start_duration = i
        elseif (start_duration != 0) && (allfemales > threshold)
            duration_of_protection = i - start_duration
            break
        end
    end

    if (start_duration != 0) && (duration_of_protection == 0)
        duration_of_protection = t - start_duration
    end


    return max_suppression, t_max_suppression, t_1e_2, t_1e_3, t_1e_4, t_1e_5,
                t_1e_6, t_1e_7, t_1e_8, t_1e_9, t_1e_10, duration_of_protection,
                transgenic_Y, transgenic_S
end



# this is the general function used for plotting run file results
# function plot_run(run::Dict; allfemales_bool = true, yle_total_bool = false,
#                     shredder_total_all_bool = false,
#                     shredder_total_male_bool = false,
#                     shredder_total_female_bool = false,
#                     Cas9_bool = false, gRNA1_bool = false,
#                     gRNA2_bool = false, shredder_bool = false,
#                     sexratio_bool = false, correlation_bool = false,
#                     y_bool = false,
#                     given_size = (700,500), given_dpi = 300)
#
#
#
#     output = plot(size = given_size, dpi = given_dpi,
#                 ylims = (-0.031,1.03))
#
#     t = length(run["genotypes"][1,:])
#     AllF = zeros(t)
#     AllM = zeros(t)
#     YLE_total = zeros(t)
#     Shredder_total_all = zeros(t)
#     Shredder_total_male = zeros(t)
#     Shredder_total_female = zeros(t)
#     Cas9 = zeros(t)
#     gRNA1 = zeros(t)
#     gRNA2 = zeros(t)
#     shredder = zeros(t)
#     sexratio = zeros(t)
#     correlation = zeros(t)
#     y = zeros(t)
#
#
#     #iterating through all generations
#     for i in 1:t
#
#         current_AllF = 0
#         current_AllM = 0
#         current_YLE_total = 0
#
#         #counting the frequency of fully functional X-shredder in all individuals
#         current_Shredder_total_all = 0
#         #counting the frequency of fully functional X-shredder in males
#         current_Shredder_total_male = 0
#         #counting the frequency of fully functional X-shredder in females
#         current_Shredder_total_female = 0
#         #counting the frequency of fully functional X-shredder in
#         #a) males only and b) with same weight for homo and heterozygous
#         current_Shredder_correlation = 0
#
#         current_both_freq = 0
#         current_Cas9 = 0
#         current_gRNA1 = 0
#         current_gRNA2 = 0
#         current_shredder = 0
#
#         current_Y = 0
#
#         #iterating through all genotypes
#         for k in 1:length(run["genotypes"][:,1])
#
#             if (genotypes_detailed[k][1] in X_chromosomes)
#                 current_AllF += run["genotypes"][k,i]
#             else
#                 current_AllM += run["genotypes"][k,i]
#             end
#
#             if ((genotypes[k][1] == "Y2") &&
#                         (genotypes[k][3] == "A2" || genotypes[k][4] == "A2" ))
#                 current_YLE_total += run["genotypes"][k,i]
#                 current_both_freq += run["genotypes"][k,i]
#             elseif (genotypes[k][1] == "Y2")
#                 current_YLE_total += run["genotypes"][k,i]
#             elseif (genotypes[k][1] == "Y1")
#                 current_Y += run["genotypes"][k,i]
#             end
#
#             if ((genotypes_detailed[k][3] == "CD") && (genotypes_detailed[k][4] == "CD"))
#                 current_Shredder_total_all += (2 * run["genotypes"][k,i])
#             elseif ((genotypes_detailed[k][3] == "CD") || (genotypes_detailed[k][4] == "CD"))
#                 current_Shredder_total_all += run["genotypes"][k,i]
#             end
#
#             if ((genotypes_detailed[k][3] == "CD") && (genotypes_detailed[k][4] == "CD") && (genotypes_detailed[k][1] ∉ X_chromosomes))
#                 current_Shredder_total_male += (2 * run["genotypes"][k,i])
#                 current_Shredder_correlation += run["genotypes"][k,i]
#             elseif (((genotypes_detailed[k][3] == "CD") || (genotypes_detailed[k][4] == "CD")) && (genotypes_detailed[k][1] ∉ X_chromosomes))
#                 current_Shredder_total_male += run["genotypes"][k,i]
#                 current_Shredder_correlation += run["genotypes"][k,i]
#             end
#
#             if ((genotypes_detailed[k][3] == "CD") && (genotypes_detailed[k][4] == "CD") && (genotypes_detailed[k][1] in X_chromosomes))
#                 current_Shredder_total_female += (2 * run["genotypes"][k,i])
#             elseif (((genotypes_detailed[k][3] == "CD") || (genotypes_detailed[k][4] == "CD")) && (genotypes_detailed[k][1] in X_chromosomes))
#                 current_Shredder_total_female += run["genotypes"][k,i]
#             end
#
#             if (occursin("A", genotypes_detailed[k][1]))
#                 current_Cas9 += run["genotypes"][k,i]
#             end
#
#             if (occursin("B", genotypes_detailed[k][1]))
#                 current_gRNA1 += run["genotypes"][k,i]
#             end
#
#             if (occursin("C", genotypes_detailed[k][3]) && occursin("C", genotypes_detailed[k][4]))
#                 current_shredder += (2 * run["genotypes"][k,i])
#             elseif (occursin("C", genotypes_detailed[k][3]) || occursin("C", genotypes_detailed[k][4]))
#                 current_shredder += run["genotypes"][k,i]
#             end
#
#             if (occursin("D", genotypes_detailed[k][3]) && occursin("D", genotypes_detailed[k][4]))
#                 current_gRNA2 += (2 * run["genotypes"][k,i])
#             elseif (occursin("D", genotypes_detailed[k][3]) || occursin("D", genotypes_detailed[k][4]))
#                 current_gRNA2 += run["genotypes"][k,i]
#             end
#         end
#
#         current_both_freq = current_both_freq / current_AllM
#         current_YLE_total = current_YLE_total / current_AllM
#         current_Shredder_correlation = current_Shredder_correlation / current_AllM
#
#         AllF[i] = current_AllF
#         YLE_total[i] = current_YLE_total
#         sexratio[i] = current_AllM / (current_AllF + current_AllM)
#         Shredder_total_all[i] = (current_Shredder_total_all /
#                                         ((current_AllM * 2) + (current_AllF*2)))
#         Shredder_total_male[i] = current_Shredder_total_male / (current_AllM * 2)
#         Shredder_total_female[i] = current_Shredder_total_female / (current_AllF * 2)
#         Cas9[i] = current_Cas9 / current_AllM
#         gRNA1[i] = current_gRNA1 / current_AllM
#         shredder[i] = current_shredder / ((current_AllM * 2) + (current_AllF * 2))
#         gRNA2[i] = current_gRNA2 / ((current_AllM * 2) + (current_AllF * 2))
#
#         y[i] = current_Y / current_AllM
#
#         correlation[i] = ((current_both_freq - (current_YLE_total * current_Shredder_correlation))/
#                        sqrt(current_YLE_total*(1 - current_YLE_total)
#                        *current_Shredder_correlation*
#                        (1 -current_Shredder_correlation)))
#     end
#
#     if allfemales_bool
#         output = plot!(1:t, AllF, lw = 4, colour=:black, linestyle=:dot,
#                     label = "Total females")
#     end
#
#     if yle_total_bool
#         output = plot!(1:t, YLE_total, lw=4, colour=:blue,
#                 label = "YLE")
#     end
#
#     if sexratio_bool
#         output = plot!(1:t, sexratio, lw=4, label = "Sex ratio (m/all)", color = 4)
#     end
#
#     if shredder_total_all_bool
#         output = plot!(1:t, Shredder_total_all, lw = 4, label = "X-shredder (all)",
#             colour = 3)
#     end
#
#     if shredder_total_male_bool
#         output = plot!(1:t, Shredder_total_male, lw = 4, label = "X-shredder (males)",
#             colour = 3)
#     end
#
#     if shredder_total_female_bool
#         output = plot!(1:t, Shredder_total_female, lw = 4, label = "X-shredder freq (females)")
#     end
#
#     if Cas9_bool
#         output = plot!(1:t, Cas9 , lw = 4, label = "Cas9")
#     end
#
#     if gRNA1_bool
#         output = plot!(1:t, gRNA1 , lw = 4, label = "gRNA1")
#     end
#
#     if shredder_bool
#         output = plot!(1:t, shredder , lw = 4, label = "shredder (single)")
#     end
#
#     if gRNA2_bool
#         output = plot!(1:t, gRNA2 , lw = 4, label = "gRNA2")
#     end
#
#     if correlation_bool
#         output = plot!(1:t, correlation, lw = 4, label = "Correlation",
#             colour = 2)
#
#             # println(correlation)
#     end
#
#     if y_bool
#         output = plot!(1:t, y, lw = 4, label = "wt Y freq")
#     end
#
#     return output
#
# end


function plot_run_dys(run::Dict; allfemales_bool = true,
                    Cas9_bool = true, gRNA1_bool = true,
                    gRNA2_bool = true, shredder_bool = true,
                    correlation_Cas9_shredder_bool = false,
                    correlation_gRNA1_shredder_bool = false,
                    given_size = (700,500), given_dpi = 300)

    output = plot(size = given_size, dpi = given_dpi)

    t = length(run["genotypes"][1,:])
    AllF = ones(t)
    AllM = ones(t)
    Cas9_dys = ones(t)
    gRNA1_dys = ones(t)
    shredder_dys = ones(t)
    gRNA2_dys = ones(t)
    correlation_Cas9_shredder = ones(t)
    correlation_gRNA1_shredder = ones(t)

    #iterating through all generations
    for i in 1:length(run["genotypes"][1,:])
        current_AllF = 0
        current_AllM = 0
        current_YLE_total = 0
        current_Cas9_dys= 0
        current_gRNA1_dys = 0
        current_shredder_dys = 0
        current_gRNA2_dys = 0

        #iterating through all genotypes
        for k in 1:length(run["genotypes"][:,1])

            if (genotypes_detailed[k][1] in X_chromosomes)
                current_AllF = current_AllF + run["genotypes"][k,i]
            else
                current_AllM = current_AllM + run["genotypes"][k,i]
            end

            if (occursin("α", genotypes_detailed[k][1]))
                current_Cas9_dys = current_Cas9_dys + run["genotypes"][k,i]
            end

            if (occursin("β", genotypes_detailed[k][1]))
                current_gRNA1_dys = current_gRNA1_dys + run["genotypes"][k,i]
            end

            if (occursin("ζ", genotypes_detailed[k][3]) && occursin("ζ", genotypes_detailed[k][4]))
                current_shredder_dys = current_shredder_dys + (2 * run["genotypes"][k,i])
            elseif (occursin("ζ", genotypes_detailed[k][3]) || occursin("ζ", genotypes_detailed[k][4]))
                current_shredder_dys = current_shredder_dys + run["genotypes"][k,i]
            end

            if (occursin("δ", genotypes_detailed[k][3]) && occursin("δ", genotypes_detailed[k][4]))
                current_gRNA2_dys = current_gRNA2_dys + (2 * run["genotypes"][k,i])
            elseif (occursin("δ", genotypes_detailed[k][3]) || occursin("δ", genotypes_detailed[k][4]))
                current_gRNA2_dys = current_gRNA2_dys + run["genotypes"][k,i]
            end

        end

    AllF[i] = current_AllF
    Cas9_dys[i] = current_Cas9_dys / current_AllM
    gRNA1_dys[i] = current_gRNA1_dys / current_AllM
    shredder_dys[i] = current_shredder_dys / ((current_AllM * 2) + (current_AllF * 2))
    gRNA2_dys[i] = current_gRNA2_dys / ((current_AllM * 2) + (current_AllF * 2))

    end

    if (correlation_Cas9_shredder_bool || correlation_gRNA1_shredder_bool)
        #iterating through all generations
        for j in 1:length(run["sperm"][1,:])

            current_dys_shredder_freq = 0
            current_dys_Cas9_freq = 0
            current_dys_gRNA1_freq = 0
            current_Cas9_shredder = 0
            current_gRNA1_shredder = 0
            current_correlation_Cas9_shredder = 0
            current_correlation_gRNA1_shredder = 0

            #iterating through all gametes
            for l in 1:length(run["sperm"][:,1])

                if (occursin("α", gametes_detailed[l][1]) && occursin("ζ", gametes_detailed[l][2]))
                    current_dys_Cas9_freq = current_dys_Cas9_freq + run["sperm"][l,j]
                    current_Cas9_shredder = current_Cas9_shredder + run["sperm"][l,j]
                elseif (occursin("α", gametes_detailed[l][1]))
                    current_dys_Cas9_freq = current_dys_Cas9_freq + run["sperm"][l,j]
                end

                if (occursin("β", gametes_detailed[l][1]) && occursin("ζ", gametes_detailed[l][2]))
                    current_dys_gRNA1_freq = current_dys_gRNA1_freq + run["sperm"][l,j]
                    current_gRNA1_shredder = current_gRNA1_shredder + run["sperm"][l,j]
                elseif (occursin("β", gametes_detailed[l][1]))
                    current_dys_gRNA1_freq = current_dys_gRNA1_freq + run["sperm"][l,j]
                end

                if (occursin("ζ", gametes_detailed[l][2]))
                    current_dys_shredder_freq = current_dys_shredder_freq + run["sperm"][l,j]
                end

            end

            # println("shredder:", current_dys_shredder_freq)
            # println("Cas9:", current_dys_Cas9_freq)
            # println("gRNA1:", current_dys_gRNA1_freq)
            # println("Cas9:", current_Cas9_shredder)
            # println("gRNA1:", current_gRNA1_shredder)

            current_correlation_Cas9_shredder = ((current_Cas9_shredder - (current_dys_Cas9_freq * current_dys_shredder_freq)) /
                                        sqrt((current_dys_Cas9_freq*(1-current_dys_Cas9_freq)*current_dys_shredder_freq*(1-current_dys_shredder_freq))))
            correlation_Cas9_shredder[j] = (current_correlation_Cas9_shredder)
            # println(current_correlation_Cas9_shredder)

            current_correlation_gRNA1_shredder = ((current_gRNA1_shredder - (current_dys_gRNA1_freq * current_dys_shredder_freq)) /
                                        sqrt((current_dys_gRNA1_freq*(1-current_dys_gRNA1_freq)*current_dys_shredder_freq*(1-current_dys_shredder_freq))))
            correlation_gRNA1_shredder[j] = (current_correlation_gRNA1_shredder)
            # println(current_correlation_gRNA1_shredder)
        end
    end


    if allfemales_bool
        output = plot!(1:t, AllF, lw = 4, label = "Total females" ,
                size = given_size, dpi = given_dpi, colour=:black, linestyle=:dot,)
    end

    if Cas9_bool
        output = plot!(1:t, Cas9_dys , lw = 4, label = "dysf Cas9",
                size = given_size, dpi = given_dpi, colour= ColorBrewer.palette("Dark2", 8)[4])
    end

    if gRNA1_bool
        output = plot!(1:t, gRNA1_dys , lw = 4, label = "dysf gRNA1",
                    size = given_size, dpi = given_dpi, colour = ColorBrewer.palette("Dark2", 8)[8])
    end

    if shredder_bool
        output = plot!(1:t, shredder_dys , lw = 4, label = "dysf shredder (single)",
                    size = given_size, dpi = given_dpi, colour = ColorBrewer.palette("Dark2", 8)[6])
    end

    if gRNA2_bool
        output = plot!(1:t, gRNA2_dys , lw = 4, label = "dysf gRNA2",
                    size = given_size, dpi = given_dpi, colour = ColorBrewer.palette("Dark2", 8)[3])
    end

    if correlation_Cas9_shredder_bool
        output = plot!(1:t, correlation_Cas9_shredder , lw = 4, label = "correl dys Cas9/shredder",size = given_size, dpi = given_dpi)
    end

    if correlation_gRNA1_shredder_bool
        output = plot!(1:t, correlation_gRNA1_shredder , lw = 4, label = "correl dys gRNA1/shredder",size = given_size, dpi = given_dpi)
    end

    return output
end


function plot_run_res(run::Dict; allfemales_bool = true, yle_total_bool = false,
                    shredder_total_all_bool = false,
                    shredder_total_male_bool = false,
                    res1_bool = true,
                    res2_bool = true,
                    res3_bool = true,
                    res_colour =:red,
                    given_size = (700,500), given_dpi = 300)

    output = plot(size = given_size, dpi = given_dpi)

    t = length(run["genotypes"][1,:])
    AllF = zeros(t)
    AllM = zeros(t)
    YLE_total = zeros(t)
    Shredder_total_all = zeros(t)
    Shredder_total_male = zeros(t)
    Res1 = zeros(t)
    Res2 = zeros(t)
    Res3 = zeros(t)



    #iterating through all generations
    for i in 1:t

        current_AllF = 0
        current_AllM = 0
        current_YLE_total = 0
        current_Shredder_total_all = 0
        current_Shredder_total_male = 0
        current_Res1 = 0
        current_Res2 = 0
        current_Res3 = 0


        #iterating through all genotypes
        for k in 1:length(run["genotypes"][:,1])

            if (genotypes_detailed[k][1] in X_chromosomes)
                current_AllF += run["genotypes"][k,i]
            else
                current_AllM += run["genotypes"][k,i]
            end

            if (genotypes[k][1] == "Y2")
                current_YLE_total += run["genotypes"][k,i]
            end

            if ((genotypes_detailed[k][3] == "CD") && (genotypes_detailed[k][4] == "CD"))
                current_Shredder_total_all += (2 * run["genotypes"][k,i])
            elseif ((genotypes_detailed[k][3] == "CD") || (genotypes_detailed[k][4] == "CD"))
                current_Shredder_total_all += run["genotypes"][k,i]
            end

            if ((genotypes_detailed[k][3] == "CD") && (genotypes_detailed[k][4] == "CD") && (genotypes_detailed[k][1] ∉ X_chromosomes))
                current_Shredder_total_male += (2 * run["genotypes"][k,i])
            elseif (((genotypes_detailed[k][3] == "CD") || (genotypes_detailed[k][4] == "CD")) && (genotypes_detailed[k][1] ∉ X_chromosomes))
                current_Shredder_total_male += run["genotypes"][k,i]
            end


            if (occursin("r1", genotypes_detailed[k][1]) && occursin("r1", genotypes_detailed[k][2]))
                current_Res1 += (2 * run["genotypes"][k,i])
            elseif (occursin("r1", genotypes_detailed[k][1]) || occursin("r1", genotypes_detailed[k][2]))
                current_Res1 += run["genotypes"][k,i]
            end

            if (occursin("r2", genotypes_detailed[k][1]) && occursin("r2", genotypes_detailed[k][2]))
                current_Res2 += (2 * run["genotypes"][k,i])
            elseif (occursin("r2", genotypes_detailed[k][1]) || occursin("r2", genotypes_detailed[k][2]))
                current_Res2 += run["genotypes"][k,i]
            end

            if (occursin("r3", genotypes_detailed[k][3]) && occursin("r3", genotypes_detailed[k][4]))
                current_Res3 += (2 * run["genotypes"][k,i])
            elseif (occursin("r3", genotypes_detailed[k][3]) || occursin("r3", genotypes_detailed[k][4]))
                current_Res3 += run["genotypes"][k,i]
            end



        end


        AllF[i] = current_AllF
        YLE_total[i] = current_YLE_total / current_AllM
        Shredder_total_all[i] = (current_Shredder_total_all /
                                        ((current_AllM * 2) + (current_AllF*2)))
        Shredder_total_male[i] = current_Shredder_total_male / (current_AllM * 2)
        Res1[i] = current_Res1 / (current_AllF * 2 + current_AllM)
        Res2[i] = current_Res2 / (current_AllF * 2 + current_AllM)
        Res3[i] = current_Res3 / (current_AllF * 2 + current_AllM * 2)

    end


    if allfemales_bool
        output = plot!(1:t, AllF, lw = 4, colour=:black, linestyle=:dot,
                    label = "Total females")
    end

    if yle_total_bool
        output = plot!(1:t, YLE_total, lw=4, colour=:blue,
                label = "YLE frequency")
    end

    if shredder_total_all_bool
        output = plot!(1:t, Shredder_total_all, lw = 4, label = "X-shredder freq (all)")
    end

    if shredder_total_male_bool
        output = plot!(1:t, Shredder_total_male, lw = 4, label = "X-shredder freq (males)")
    end

    if res1_bool
        output = plot!(1:t, Res1, lw = 4, label = "editing resist freq", colour = res_colour)
    end

    if res2_bool
        output = plot!(1:t, Res2, lw = 4, label = "shredding resist freq", colour = res_colour)
    end

    if res3_bool
        output = plot!(1:t, Res3, lw = 4, label = "homing resist freq", colour = res_colour)
    end

    return output

end




function plot_run_comprehensive(run::Dict; allfemales_bool = true, yle_total_bool = false,
                    shredder_total_all_bool = false,
                    shredder_total_male_bool = false,
                    shredder_total_female_bool = false,
                    Cas9_bool = false, gRNA1_bool = false,
                    gRNA2_bool = false, shredder_bool = false,
                    dCas9_bool = false, dgRNA1_bool = false,
                    dgRNA2_bool = false, dshredder_bool = false,
                    res1_bool = false,
                    res2_bool = false,
                    res3_bool = false,
                    sexratio_bool = false, correlation_bool = false,
                    y_bool = false,
                    load_easy_bool = false,
                    load_demanding_bool = false,
                    inverse_load_bool = false,
                    given_size = (700,500), given_dpi = 300,
                    log_safe = false,
                    allfemales_log = false,
                    second_y_log_lims = (1e-8, 1.5),
                    second_y_ticks = [1, 1e-2, 1e-4, 1e-6, 1e-8],
                    dCas9_colour = colour_dCas9,
                    dgRNA1_colour = colour_dgRNA1 ,
                    dshredder_colour =  colour_dShredder,
                    dgRNA2_colour = colour_dgRNA2,
                    load_colour = colour_load,
                    initial_size = 1.0,
                    print_load = false,
                    compare_load_calcuations = false,
                    load_zygotes_bool = false,
                    load_zygotes_bool2 = false)

    global current_Parameters, current_matrices

    output = plot(size = given_size, dpi = given_dpi,
                ylims = (-0.031,1.03))

    t = length(run["genotypes"][1,:])
    AllF = zeros(t)
    AllM = zeros(t)
    YLE_total = zeros(t)
    Shredder_total_all = zeros(t)
    Shredder_total_male = zeros(t)
    Shredder_total_female = zeros(t)
    Cas9 = zeros(t)
    gRNA1 = zeros(t)
    gRNA2 = zeros(t)
    Cas9_dys = ones(t)
    gRNA1_dys = ones(t)
    shredder_dys = ones(t)
    gRNA2_dys = ones(t)
    shredder = zeros(t)
    sexratio = zeros(t)
    correlation = zeros(t)
    Res1 = zeros(t)
    Res2 = zeros(t)
    Res3 = zeros(t)
    y = zeros(t)
    load_easy = ones(t)
    load_demanding = ones(t)


    Rm = current_Parameters["Rm"]
    theta = current_Parameters["theta"]

    #number of eggs ("f")
    f = (Rm * 2.0) / theta

    load_zygotes = zeros(t)
    load_zygotes2 = zeros(t)

    #alpha is a constant determining the density dependent mortality (half maximal)
    alpha = initial_size * f / (Rm - 1.0)


    #iterating through all generations
    for i in 1:t

        current_AllF = 0
        current_AllM = 0
        current_YLE_total = 0

        current_AllF_zygotes = 0
        current_AllF_x_selection_zygotes = 0

        #counting the frequency of fully functional X-shredder in all individuals
        current_Shredder_total_all = 0
        #counting the frequency of fully functional X-shredder in males
        current_Shredder_total_male = 0
        #counting the frequency of fully functional X-shredder in females
        current_Shredder_total_female = 0
        #counting the frequency of fully functional X-shredder in
        #a) males only and b) with same weight for homo and heterozygous
        current_Shredder_correlation = 0

        current_both_freq = 0
        current_Cas9 = 0
        current_gRNA1 = 0
        current_gRNA2 = 0
        current_shredder = 0

        current_Cas9_dys= 0
        current_gRNA1_dys = 0
        current_shredder_dys = 0
        current_gRNA2_dys = 0

        current_Res1 = 0
        current_Res2 = 0
        current_Res3 = 0

        current_Y = 0

        #counting the number of females in alternative scenario without
        #editing, homing etc
        wildtype_AllF = 0

        #iterating through all genotypes
        for k in 1:length(run["genotypes"][:,1])

            if (genotypes_detailed[k][1] in X_chromosomes)
                current_AllF += run["genotypes"][k,i]
                wildtype_AllF += run["wildtype"][k,i]
                current_AllF_zygotes += run["zygotes"][k,i]
                current_AllF_x_selection_zygotes += run["zygotes"][k,i] * current_matrices["selection_matrix"][k]
            else
                current_AllM += run["genotypes"][k,i]
            end


            if ((genotypes[k][1] == "Y2") &&
                        (genotypes[k][3] == "A2" || genotypes[k][4] == "A2" ))
                current_YLE_total += run["genotypes"][k,i]
                current_both_freq += run["genotypes"][k,i]
            elseif (genotypes[k][1] == "Y2")
                current_YLE_total += run["genotypes"][k,i]
            elseif (genotypes[k][1] == "Y1")
                current_Y += run["genotypes"][k,i]
            end

            if ((genotypes_detailed[k][3] == "CD") && (genotypes_detailed[k][4] == "CD"))
                current_Shredder_total_all += (2 * run["genotypes"][k,i])
            elseif ((genotypes_detailed[k][3] == "CD") || (genotypes_detailed[k][4] == "CD"))
                current_Shredder_total_all += run["genotypes"][k,i]
            end

            if ((genotypes_detailed[k][3] == "CD") && (genotypes_detailed[k][4] == "CD") && (genotypes_detailed[k][1] ∉ X_chromosomes))
                current_Shredder_total_male += (2 * run["genotypes"][k,i])
                current_Shredder_correlation += run["genotypes"][k,i]
            elseif (((genotypes_detailed[k][3] == "CD") || (genotypes_detailed[k][4] == "CD")) && (genotypes_detailed[k][1] ∉ X_chromosomes))
                current_Shredder_total_male += run["genotypes"][k,i]
                current_Shredder_correlation += run["genotypes"][k,i]
            end

            if ((genotypes_detailed[k][3] == "CD") && (genotypes_detailed[k][4] == "CD") && (genotypes_detailed[k][1] in X_chromosomes))
                current_Shredder_total_female += (2 * run["genotypes"][k,i])
            elseif (((genotypes_detailed[k][3] == "CD") || (genotypes_detailed[k][4] == "CD")) && (genotypes_detailed[k][1] in X_chromosomes))
                current_Shredder_total_female += run["genotypes"][k,i]
            end

            if (occursin("A", genotypes_detailed[k][1]))
                current_Cas9 += run["genotypes"][k,i]
            end

            if (occursin("B", genotypes_detailed[k][1]))
                current_gRNA1 += run["genotypes"][k,i]
            end

            if (occursin("C", genotypes_detailed[k][3]) && occursin("C", genotypes_detailed[k][4]))
                current_shredder += (2 * run["genotypes"][k,i])
            elseif (occursin("C", genotypes_detailed[k][3]) || occursin("C", genotypes_detailed[k][4]))
                current_shredder += run["genotypes"][k,i]
            end

            if (occursin("D", genotypes_detailed[k][3]) && occursin("D", genotypes_detailed[k][4]))
                current_gRNA2 += (2 * run["genotypes"][k,i])
            elseif (occursin("D", genotypes_detailed[k][3]) || occursin("D", genotypes_detailed[k][4]))
                current_gRNA2 += run["genotypes"][k,i]
            end


            if (occursin("α", genotypes_detailed[k][1]))
                current_Cas9_dys = current_Cas9_dys + run["genotypes"][k,i]
            end

            if (occursin("β", genotypes_detailed[k][1]))
                current_gRNA1_dys = current_gRNA1_dys + run["genotypes"][k,i]
            end

            if (occursin("ζ", genotypes_detailed[k][3]) && occursin("ζ", genotypes_detailed[k][4]))
                current_shredder_dys = current_shredder_dys + (2 * run["genotypes"][k,i])
            elseif (occursin("ζ", genotypes_detailed[k][3]) || occursin("ζ", genotypes_detailed[k][4]))
                current_shredder_dys = current_shredder_dys + run["genotypes"][k,i]
            end

            if (occursin("δ", genotypes_detailed[k][3]) && occursin("δ", genotypes_detailed[k][4]))
                current_gRNA2_dys = current_gRNA2_dys + (2 * run["genotypes"][k,i])
            elseif (occursin("δ", genotypes_detailed[k][3]) || occursin("δ", genotypes_detailed[k][4]))
                current_gRNA2_dys = current_gRNA2_dys + run["genotypes"][k,i]
            end

            if (occursin("r1", genotypes_detailed[k][1]) && occursin("r1", genotypes_detailed[k][2]))
                current_Res1 += (2 * run["genotypes"][k,i])
            elseif (occursin("r1", genotypes_detailed[k][1]) || occursin("r1", genotypes_detailed[k][2]))
                current_Res1 += run["genotypes"][k,i]
            end

            if (occursin("r2", genotypes_detailed[k][1]) && occursin("r2", genotypes_detailed[k][2]))
                current_Res2 += (2 * run["genotypes"][k,i])
            elseif (occursin("r2", genotypes_detailed[k][1]) || occursin("r2", genotypes_detailed[k][2]))
                current_Res2 += run["genotypes"][k,i]
            end

            if (occursin("r3", genotypes_detailed[k][3]) && occursin("r3", genotypes_detailed[k][4]))
                current_Res3 += (2 * run["genotypes"][k,i])
            elseif (occursin("r3", genotypes_detailed[k][3]) || occursin("r3", genotypes_detailed[k][4]))
                current_Res3 += run["genotypes"][k,i]
            end


        end

        current_both_freq = current_both_freq / current_AllM
        current_YLE_total = current_YLE_total / current_AllM
        current_Shredder_correlation = current_Shredder_correlation / current_AllM

        AllF[i] = current_AllF
        YLE_total[i] = current_YLE_total
        sexratio[i] = current_AllF / (current_AllF + current_AllM)
        Shredder_total_all[i] = (current_Shredder_total_all /
                                        ((current_AllM * 2) + (current_AllF*2)))
        Shredder_total_male[i] = current_Shredder_total_male / (current_AllM * 2)
        Shredder_total_female[i] = current_Shredder_total_female / (current_AllF * 2)
        Cas9[i] = current_Cas9 / current_AllM
        gRNA1[i] = current_gRNA1 / current_AllM
        shredder[i] = current_shredder / ((current_AllM * 2) + (current_AllF * 2))
        gRNA2[i] = current_gRNA2 / ((current_AllM * 2) + (current_AllF * 2))


        Cas9_dys[i] = current_Cas9_dys / current_AllM
        gRNA1_dys[i] = current_gRNA1_dys / current_AllM
        shredder_dys[i] = current_shredder_dys / ((current_AllM * 2) + (current_AllF * 2))
        gRNA2_dys[i] = current_gRNA2_dys / ((current_AllM * 2) + (current_AllF * 2))

        if (load_easy_bool && i<t)
            load_easy[i+1] = 0.5 * current_AllF * f * theta * alpha / (alpha + (current_AllF * f))

            load_easy[i] = (current_AllF) / load_easy[i]
            println("Easy load at ", i, " is ", 1 - load_easy[i])

            if (print_load)
                println(load_easy[i])
            end

        else
            load_easy[i] = NaN
        end




        load_demanding[i] = current_AllF / wildtype_AllF

        if (print_load)
            println(1 - load_demanding[i])
        end

        if (compare_load_calcuations)
            if (load_demanding[i] != load_easy[i])
                println("Warning: load calculation methods differ for")
                println(i)
                println(load_easy[i])
                println(load_demanding[i])
            end
        end


        if (load_zygotes_bool  && i<t)
            #println(i)
            #load_zygotes[i+1] = 0.5 * current_AllF * f * theta * alpha / (alpha + (current_AllF * f))
            #load_zygotes[i + 1] = 0.5 * current_AllF * f
            load_zygotes[i + 1] = (current_AllF_zygotes * theta * alpha / (alpha + (sum(
                        run["zygotes"][:,i])))) *
                                        f * 0.5


            #println("has allF_zygotes:", current_AllF_zygotes)
            #println("has current_AllF_x_selection_zygotes:", current_AllF_x_selection_zygotes)

            #println("Next generation has", load_zygotes[i+1], "zygotes")
            load_zygotes[i] = 1 - ((current_AllF_x_selection_zygotes) / load_zygotes[i])
            #println("Zygotes load 1 at ", i, " is ", load_zygotes[i])
            #println("This generation has load of", load_zygotes[i])
        elseif (load_zygotes_bool)
            load_zygotes[i] = 1 - ((current_AllF_x_selection_zygotes) / load_zygotes[i])
        end

        if (load_zygotes_bool2)

            load_zygotes2[i] = 1 - 2 * (current_AllF_zygotes/sum(
                        run["zygotes"][:,i])) * ((current_AllF_x_selection_zygotes) / current_AllF_zygotes)

            #println("Zygotes load2 at ", i, " is ", load_zygotes2[i])

        end


        if log_safe
            if Cas9_dys[i] == 0
                Cas9_dys[i] = NaN
            end
            if gRNA1_dys[i] == 0
                gRNA1_dys[i] = NaN
            end
            if shredder_dys[i] == 0
                shredder_dys[i] = NaN
            end
            if gRNA2_dys[i] == 0
                gRNA2_dys[i] = NaN
            end
        end

        Res1[i] = current_Res1 / (current_AllF * 2 + current_AllM)
        Res2[i] = current_Res2 / (current_AllF * 2 + current_AllM)
        Res3[i] = current_Res3 / (current_AllF * 2 + current_AllM * 2)


        y[i] = current_Y / current_AllM



        correlation[i] = ((current_both_freq - (current_YLE_total * current_Shredder_correlation))/
                       sqrt(current_YLE_total*(1 - current_YLE_total)
                       *current_Shredder_correlation*
                       (1 - current_Shredder_correlation)))
    end


    if sexratio_bool
        output = plot!(1:t, sexratio, lw=4, label = "Sex ratio (m/all)", color = colour_sexratio)
    end

    if shredder_total_female_bool
        output = plot!(1:t, Shredder_total_female, lw = 4, label = "X-shredder freq (females)")
    end

    if Cas9_bool
        output = plot!(1:t, Cas9 , lw = 4, label = "Cas9")
    end

    if gRNA1_bool
        output = plot!(1:t, gRNA1 , lw = 4, label = "gRNA1")
    end

    if shredder_bool
        output = plot!(1:t, shredder , lw = 4, label = "shredder (single)")
    end

    if gRNA2_bool
        output = plot!(1:t, gRNA2 , lw = 4, label = "gRNA2")
    end


    if dgRNA1_bool
        output = plot!(1:t, gRNA1_dys , lw = 4, label = "dysf gRNA1",
                    size = given_size, dpi = given_dpi, colour = dgRNA1_colour,
                    linestyle =:dot)
    end

    if dshredder_bool
        output = plot!(1:t, shredder_dys , lw = 4, label = "dysf shredder (single)",
                    size = given_size, dpi = given_dpi, colour = dshredder_colour,
                    linestyle =:dash)
    end

    if dgRNA2_bool
        output = plot!(1:t, gRNA2_dys , lw = 4, label = "dysf gRNA2",
                    size = given_size, dpi = given_dpi, colour = dgRNA2_colour,
                    linestyle =:dot)
    end

    if y_bool
        output = plot!(1:t, y, lw = 4, label = "wt Y freq")
    end

    if yle_total_bool
        output = plot!(1:t, YLE_total, lw=4, colour = colour_YLE,
                label = "YLE")
    end

    if shredder_total_all_bool
        output = plot!(1:t, Shredder_total_all, lw = 4, label = "X-shredder (all)",
            colour = colour_XShredder)
    end

    if shredder_total_male_bool
        output = plot!(1:t, Shredder_total_male, lw = 4, label = "X-shredder (males)",
            colour = colour_XShredder)
    end


    if correlation_bool
        output = plot!(1:t, correlation, lw = 4, label = "Correlation",
            colour = colour_correlation)
            # println(correlation)
    end

    if load_easy_bool
        output = plot!(1:t, 1 .- load_easy, lw = 4, label = "load_easy", colour = load_colour,
        linestyle =:solid)
    end

    if load_demanding_bool
        output = plot!(1:t, load_demanding, lw = 4, label = "load_demanding", colour = load_colour,
        linestyle =:solid)
    end

    if inverse_load_bool
        output = plot!(1:t, 1 .- load_demanding, lw = 4, label = "inverse load", colour = load_colour,
        linestyle =:solid)
    end

    if load_zygotes_bool
        output = plot!(1:t, load_zygotes, lw = 4, label = "load zygotes", colour = load_colour,
        linestyle =:solid)
    end

    if load_zygotes_bool2
        output = plot!(1:t, load_zygotes2, lw = 4, label = "load zygotes2", colour = load_colour,
        linestyle =:solid)
    end

    if dCas9_bool
        output = plot!(1:t, Cas9_dys , lw = 4, label = "dysf Cas9",
                size = given_size, dpi = given_dpi, colour= dCas9_colour,
                linestyle =:dash)
    end

    if res1_bool
        output = plot!(1:t, Res1, lw = 4, label = "editing resist freq", colour = colour_resistances,
        linestyle =:solid)
    end

    if res2_bool
        output = plot!(1:t, Res2, lw = 4, label = "shredding resist freq", colour = colour_resistances,
        linestyle =:solid)
    end

    if res3_bool
        output = plot!(1:t, Res3, lw = 4, label = "homing resist freq", colour = colour_resistances,
        linestyle =:solid)
    end

    if allfemales_bool && !allfemales_log
        output = plot!(1:t, AllF, lw = 4, colour = colour_allfemales, linestyle=:dot,
                    label = "Total females")
    end

    if allfemales_bool && allfemales_log
        output = plot!(twinx(), 1:t, AllF, lw = 4, colour = colour_allfemales, linestyle=:dot,
                    label = "Total females", yaxis=:log, ylims = second_y_log_lims,
                    yticks = second_y_ticks, grid = false)
    end

    output = plot!(xlims = (0,t), grid = false)

    return output
end



function plot_resistance_scan(scan, lower_lim, upper_lim; second_lower_lim = 0.0)

    sup_column = findfirst(x->x=="max_sup", scan[1,:])
    Y_column = findfirst(x->x=="Y", scan[1,:])
    S_column = findfirst(x->x=="S", scan[1,:])


    output = plot(xlims=(-0.02,1.02), ylims=(-0.02,1.02))
    plot!(Shape([0, lower_lim, lower_lim, 0], [-1, -1, 2, 2]), colour=:grey, opacity= 0.3)

    if second_lower_lim != 0.0
        plot!(Shape([0, second_lower_lim, second_lower_lim, 0], [-1, -1, 2, 2]), colour=:grey, opacity= 0.3)
    end

    plot!(Shape([upper_lim, 1, 1, upper_lim], [-1, -1, 2, 2]), colour=:grey, opacity= 0.3)
    plot!(scan[2:end,1], scan[2:end, sup_column], lw = 4, label = "min population size", color =:black, linestyle =:dot,
    dpi = 300, size = (325,150))
    plot!(scan[2:end,1], scan[2:end, Y_column], lw = 4, label = "transgenic Y", colour = colour_YLE)
    plot!(scan[2:end,1], scan[2:end, S_column], lw = 4, label = "transgenic autosomes", colour = colour_XShredder)
    plot!(legend=false, grid=false)

    return output
end

# this function takes a given set of parameter values and plugs them into
# the different tables to convert those symbolic terms into numeric values
function apply_parameters_set(given_Parameters_set)

    global current_Parameters, current_matrices, selection_matrix, sperm_matrix,
            egg_matrix, editing_matrix, homing_matrix, mutation_matrix

    current_Parameters = deepcopy(given_Parameters_set)

    # if (current_Parameters["e_e"] + current_Parameters["er_1"]) > 1

    #     println("Your parameters are inconsistent!")
    #     println("e_e + er_1 > 1")
    # end
    #
    # if (current_Parameters["e_s"] + current_Parameters["er_2"]) > 1

    #     println("Your parameters are inconsistent!")
    #     println("e_s + er_2 > 1")
    # end
    #
    # if (current_Parameters["e_h"] + current_Parameters["er_3"]) > 1

    #     println("Your parameters are inconsistent!")
    #     println("e_h + er_3 > 1")
    # end

    selection_matrix0 = BigFloat.(float.(selection_matrix.evalf(subs = current_Parameters, given_precision)))

    sperm_matrix0 = deepcopy(sperm_matrix_num)

    for i in 1:length(sperm_matrix[:,1])
        for j in sperm_symbols[i]
            sperm_matrix0[i,j] = BigFloat(float.(sperm_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end


    egg_matrix0 = deepcopy(egg_matrix_num)

    for i in 1:length(egg_matrix[:,1])
        for j in egg_symbols[i]
            egg_matrix0[i,j] = BigFloat(float.(egg_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end


    editing_matrix0 = deepcopy(editing_matrix_num)

    for i in 1:length(editing_matrix[:,1])
        for j in editing_symbols[i]
            editing_matrix0[i,j] = BigFloat(float.(editing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end



    homing_matrix0 = deepcopy(homing_matrix_num)

    for i in 1:length(homing_matrix[:,1])
        for j in homing_symbols[i]
            homing_matrix0[i,j] = BigFloat(float.(homing_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end


    mutation_matrix0 = deepcopy(mutation_matrix_num)

    for i in 1:length(mutation_matrix[:,1])
        for j in mutation_symbols[i]
            mutation_matrix0[i,j] = BigFloat(float.(mutation_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end


    recombination_matrix0 = deepcopy(recombination_matrix_num)

    for i in 1:length(recombination_matrix[:,1])
        for j in recombination_symbols[i]
            recombination_matrix0[i,j] = BigFloat(float.(recombination_matrix0[i,j].evalf(subs = current_Parameters, given_precision)))
        end
    end


    current_matrices = Dict(
        "selection_matrix" => selection_matrix0,
        "egg_matrix" => egg_matrix0,
        "sperm_matrix" => sperm_matrix0,
        "editing_matrix" => editing_matrix0,
        "homing_matrix" => homing_matrix0,
        "mutation_matrix" => mutation_matrix0,
        "recombination_matrix" => recombination_matrix0,
    )

end


# this function generates a population in equilibrium with resistance
# frequencies as they are passed
function generate_population(;r1 = 0.0, r2 = 0.0, r3 = 0.0)

    output = BigFloat.(zeros(length(genotypes)))
    e = (1.0 - r1)
    f = (1.0 - r2)

    Y = ["Y1" (1.0)]
    X = ["X1" (e*f);
         "X2" (e*r2);
         "X5" (r1*f);
         "X6" (r1*r2);]
    A = ["A1" (1.0 - r3);
         "A6" (r3);]

    for i in 1:length(Y[:,1])
        for j in 1:length(X[:,1])
            for k in 1:length(A[:,1])
                for l in 1:length(A[:,1])
                    i_genotype = return_i([Y[i,1] X[j,1] A[k,1] A[l,1]], genotypes)
                    output[i_genotype] += BigFloat(Y[i,2] * X[j,2] * A[k,2] * A[l,2])
                    # println(i_genotype)
                    # println((Y[i,2] * X[j,2] * A[k,2] * A[l,2]))
                end
            end
        end
    end

    for i in 1:length(X[:,1])
        for j in 1:length(X[:,1])
            for k in 1:length(A[:,1])
                for l in 1:length(A[:,1])
                    i_genotype = return_i([X[i,1] X[j,1] A[k,1] A[l,1]], genotypes)
                    output[i_genotype] += BigFloat(X[i,2] * X[j,2] * A[k,2] * A[l,2])
                    # println(i_genotype)
                    # println((Y[i,2] * X[j,2] * A[k,2] * A[l,2]))
                end
            end
        end
    end

    return output

end

function duo_timecourse(t, genotype_vector1, genotype_vector2;
                            rel_size = 1.0, migration_rate = 0.001,
                            initial_size = 1.0)
    global current_Parameters, current_matrices

    store_genotypes1 = deepcopy(genotype_vector1)
    # store_zygotes1 = zeros(length(genotype_vector1))
    # store_eggs1 = zeros(length(gametes_detailed))
    # store_sperm1 = zeros(length(gametes_detailed))

    store_genotypes2 = deepcopy(genotype_vector2)
    # store_zygotes2 = zeros(length(genotype_vector2))
    # store_eggs2 = zeros(length(gametes_detailed))
    # store_sperm2 = zeros(length(gametes_detailed))

    Rm = current_Parameters["Rm"]
    theta = current_Parameters["theta"]

    #number of eggs ("f")
    f = (Rm * 2.0) / theta

    #alpha is a constant determining the density dependent mortality (half maximal)
    alpha = initial_size * f / (Rm - 1.0)

    for i in 1:t

        post_migration1 = zeros(length(genotype_vector1))

        for j in 1:length(post_migration1)
            post_migration1[j] = ((genotype_vector1[j] * (1 - migration_rate))
                                + (genotype_vector2[j] * migration_rate))
        end

        post_migration2 = zeros(length(genotype_vector2))

        for j in 1:length(post_migration2)
            post_migration2[j] = ((genotype_vector2[j] * (1 - migration_rate))
                                + (genotype_vector1[j] * migration_rate))
        end

        genotype_vector1 = post_migration1
        genotype_vector2 = post_migration2

        #calculating prcesses in population 1
        genotype_vector1 = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector1, mutation_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector1, homing_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector1, editing_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector1, recombination_matrix_helper)

        sperm1 = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector1, sperm_matrix_helper)
        total_sperm1 = sum(sperm1)
        sperm1 = sperm1./total_sperm1

        eggs1 = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector1, egg_matrix_helper)
        eggs1 = eggs1 .* f

        # store_sperm1 = hcat(store_sperm1, sperm1)
        # store_eggs1 = hcat(store_eggs1, eggs1)

        genotype_vector1 = create_zygote_vector(sperm1, eggs1)

        # store_zygotes1 = hcat(store_zygotes1, genotype_vector1)

        #println("RUN ", i,)
        #println(sperm)
        #println(eggs)

        total_zygotes1 = sum(genotype_vector1)

        genotype_vector1 = genotype_vector1 .* (theta * alpha / (alpha + total_zygotes1))

        genotype_vector1 = current_matrices["selection_matrix"] .* genotype_vector1

        #calculating processes in population 2
        genotype_vector2 = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector2, mutation_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector2, homing_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector2, editing_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector2, recombination_matrix_helper)

        sperm2 = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector2, sperm_matrix_helper)
        total_sperm2 = sum(sperm2)
        sperm2 = sperm2./total_sperm2

        eggs2 = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector2, egg_matrix_helper)
        eggs2 = eggs2 .* f

        # store_sperm2 = hcat(store_sperm2, sperm2)
        # store_eggs2 = hcat(store_eggs2, eggs2)

        genotype_vector2 = create_zygote_vector(sperm2, eggs2)

        # store_zygotes2 = hcat(store_zygotes2, genotype_vector2)

        #println("RUN ", i,)
        #println(sperm)
        #println(eggs)

        total_zygotes2 = sum(genotype_vector2)

        genotype_vector2 = genotype_vector2 .* (theta * alpha / (alpha + total_zygotes2))

        genotype_vector2 = current_matrices["selection_matrix"] .* genotype_vector2


        store_genotypes1 = hcat(store_genotypes1, genotype_vector1)

        store_genotypes2 = hcat(store_genotypes2, genotype_vector2)

    end

    run1 = Dict(
        "genotypes" => store_genotypes1,
        # "zygotes" => store_zygotes1,
        # "eggs" => store_eggs1,
        # "sperm" => store_sperm1,
    )

    run2 = Dict(
        "genotypes" => store_genotypes2,
        # "zygotes" => store_zygotes2,
        # "eggs" => store_eggs2,
        # "sperm" => store_sperm2,
    )

    return run1, run2

end


function duo_timecourse_uni(t, genotype_vector1, genotype_vector2;
                            rel_size = 1.0, migration_rate = 0.001,
                            initial_size = 1.0)
    global current_Parameters, current_matrices

    store_genotypes1 = deepcopy(genotype_vector1)
    # store_zygotes1 = zeros(length(genotype_vector1))
    # store_eggs1 = zeros(length(gametes_detailed))
    # store_sperm1 = zeros(length(gametes_detailed))

    store_genotypes2 = deepcopy(genotype_vector2)
    # store_zygotes2 = zeros(length(genotype_vector2))
    # store_eggs2 = zeros(length(gametes_detailed))
    # store_sperm2 = zeros(length(gametes_detailed))

    Rm = current_Parameters["Rm"]
    theta = current_Parameters["theta"]

    #number of eggs ("f")
    f = (Rm * 2.0) / theta

    #alpha is a constant determining the density dependent mortality (half maximal)
    alpha = initial_size * f / (Rm - 1.0)

    for i in 1:t

        post_migration1 = zeros(length(genotype_vector1))

        for j in 1:length(post_migration1)
            post_migration1[j] = ((genotype_vector1[j] * (1 - migration_rate)))
        end

        post_migration2 = zeros(length(genotype_vector2))

        for j in 1:length(post_migration2)
            post_migration2[j] = ((genotype_vector2[j])
                                + (genotype_vector1[j] * migration_rate))
        end

        # genotype_vector1 = post_migration1
        genotype_vector2 = post_migration2

        #calculating prcesses in population 1
        genotype_vector1 = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector1, mutation_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector1, homing_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector1, editing_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector1, recombination_matrix_helper)

        sperm1 = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector1, sperm_matrix_helper)
        total_sperm1 = sum(sperm1)
        sperm1 = sperm1./total_sperm1

        eggs1 = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector1, egg_matrix_helper)
        eggs1 = eggs1 .* f

        # store_sperm1 = hcat(store_sperm1, sperm1)
        # store_eggs1 = hcat(store_eggs1, eggs1)

        genotype_vector1 = create_zygote_vector(sperm1, eggs1)

        # store_zygotes1 = hcat(store_zygotes1, genotype_vector1)

        #println("RUN ", i,)
        #println(sperm)
        #println(eggs)

        total_zygotes1 = sum(genotype_vector1)

        genotype_vector1 = genotype_vector1 .* (theta * alpha / (alpha + total_zygotes1))

        genotype_vector1 = current_matrices["selection_matrix"] .* genotype_vector1

        #calculating processes in population 2
        genotype_vector2 = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector2, mutation_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector2, homing_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector2, editing_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector2, recombination_matrix_helper)

        sperm2 = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector2, sperm_matrix_helper)
        total_sperm2 = sum(sperm2)
        sperm2 = sperm2./total_sperm2

        eggs2 = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector2, egg_matrix_helper)
        eggs2 = eggs2 .* f

        # store_sperm2 = hcat(store_sperm2, sperm2)
        # store_eggs2 = hcat(store_eggs2, eggs2)

        genotype_vector2 = create_zygote_vector(sperm2, eggs2)

        # store_zygotes2 = hcat(store_zygotes2, genotype_vector2)

        #println("RUN ", i,)
        #println(sperm)
        #println(eggs)

        total_zygotes2 = sum(genotype_vector2)

        genotype_vector2 = genotype_vector2 .* (theta * alpha / (alpha + total_zygotes2))

        genotype_vector2 = current_matrices["selection_matrix"] .* genotype_vector2


        store_genotypes1 = hcat(store_genotypes1, genotype_vector1)

        store_genotypes2 = hcat(store_genotypes2, genotype_vector2)

    end

    run1 = Dict(
        "genotypes" => store_genotypes1,
        # "zygotes" => store_zygotes1,
        # "eggs" => store_eggs1,
        # "sperm" => store_sperm1,
    )

    run2 = Dict(
        "genotypes" => store_genotypes2,
        # "zygotes" => store_zygotes2,
        # "eggs" => store_eggs2,
        # "sperm" => store_sperm2,
    )

    return run1, run2

end



#this function is the same as the previous one but returns both max suppression and
#duration of protection
#and also already uses the optimised multiplication function
function duo_timecourse_parameterscan2(t, genotype_vector1, genotype_vector2;
                            rel_size = 1.0, migration_rate = 0.001,
                            initial_size = 1.0, threshold = 1e-2)
    global current_Parameters, current_matrices

    max_suppression = 1
    duration_of_protection = 0
    start_duration = 0


    sex_breakpoint = 0
    for m in 1:length(genotypes_detailed)
        if genotypes_detailed[m][1] in X_chromosomes
            sex_breakpoint = m
            break
        end
    end

    # store_genotypes1 = deepcopy(genotype_vector1)
    # store_zygotes1 = zeros(length(genotype_vector1))
    # store_eggs1 = zeros(length(gametes_detailed))
    # store_sperm1 = zeros(length(gametes_detailed))

    # store_genotypes2 = deepcopy(genotype_vector2)
    # store_zygotes2 = zeros(length(genotype_vector2))
    # store_eggs2 = zeros(length(gametes_detailed))
    # store_sperm2 = zeros(length(gametes_detailed))

    Rm = current_Parameters["Rm"]
    theta = current_Parameters["theta"]

    #number of eggs ("f")
    f = (Rm * 2.0) / theta

    #alpha is a constant determining the density dependent mortality (half maximal)
    alpha = initial_size * f / (Rm - 1.0)

    for i in 1:t

        post_migration1 = zeros(length(genotype_vector1))

        for j in 1:length(post_migration1)
            post_migration1[j] = ((genotype_vector1[j] * (1 - migration_rate))
                                + (genotype_vector2[j] * migration_rate))
        end

        post_migration2 = zeros(length(genotype_vector2))

        for j in 1:length(post_migration2)
            post_migration2[j] = ((genotype_vector2[j] * (1 - migration_rate))
                                + (genotype_vector1[j] * migration_rate))
        end

        genotype_vector1 = post_migration1
        genotype_vector2 = post_migration2

        #calculating prcesses in population 1
        genotype_vector1 = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector1, mutation_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector1, homing_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector1, editing_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector1, recombination_matrix_helper)

        sperm1 = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector1, sperm_matrix_helper)
        total_sperm1 = sum(sperm1)
        sperm1 = sperm1./total_sperm1

        eggs1 = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector1, egg_matrix_helper)
        eggs1 = eggs1 .* f

        # store_sperm1 = hcat(store_sperm1, sperm1)
        # store_eggs1 = hcat(store_eggs1, eggs1)

        genotype_vector1 = create_zygote_vector(sperm1, eggs1)

        # store_zygotes1 = hcat(store_zygotes1, genotype_vector1)

        #println("RUN ", i,)
        #println(sperm)
        #println(eggs)

        total_zygotes1 = sum(genotype_vector1)

        genotype_vector1 = genotype_vector1 .* (theta * alpha / (alpha + total_zygotes1))

        genotype_vector1 = current_matrices["selection_matrix"] .* genotype_vector1

        #calculating processes in population 2
        genotype_vector2 = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector2, mutation_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector2, homing_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector2, editing_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector2, recombination_matrix_helper)

        sperm2 = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector2, sperm_matrix_helper)
        total_sperm2 = sum(sperm2)
        sperm2 = sperm2./total_sperm2

        eggs2 = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector2, egg_matrix_helper)
        eggs2 = eggs2 .* f

        # store_sperm2 = hcat(store_sperm2, sperm2)
        # store_eggs2 = hcat(store_eggs2, eggs2)

        genotype_vector2 = create_zygote_vector(sperm2, eggs2)

        # store_zygotes2 = hcat(store_zygotes2, genotype_vector2)

        #println("RUN ", i,)
        #println(sperm)
        #println(eggs)

        total_zygotes2 = sum(genotype_vector2)

        genotype_vector2 = genotype_vector2 .* (theta * alpha / (alpha + total_zygotes2))

        genotype_vector2 = current_matrices["selection_matrix"] .* genotype_vector2


        # store_genotypes1 = hcat(store_genotypes1, genotype_vector1)
        #
        # store_genotypes2 = hcat(store_genotypes2, genotype_vector2)

        allfemales = 0
        for j in sex_breakpoint:length(genotypes_detailed)
            allfemales = allfemales + genotype_vector1[j]
        end

        if allfemales < max_suppression
            max_suppression = allfemales
        end


        if (start_duration == 0) && (allfemales < threshold)
            start_duration = i
        elseif (start_duration != 0) && (allfemales > threshold)
            duration_of_protection = i - start_duration
            break
        end

    end

    # run1 = Dict(
    #     "genotypes" => store_genotypes1,
    #     # "zygotes" => store_zygotes1,
    #     # "eggs" => store_eggs1,
    #     # "sperm" => store_sperm1,
    # )
    #
    # run2 = Dict(
    #     "genotypes" => store_genotypes2,
    #     # "zygotes" => store_zygotes2,
    #     # "eggs" => store_eggs2,
    #     # "sperm" => store_sperm2,
    # )

    if (start_duration != 0) && (duration_of_protection == 0)
        duration_of_protection = t - start_duration
    end

    return max_suppression, duration_of_protection

end

function duo_timecourse_parameterscan2_unidirectional_target(t, genotype_vector1, genotype_vector2;
                            rel_size = 1.0, migration_rate = 0.001,
                            initial_size = 1.0, threshold = 1e-2)
    global current_Parameters, current_matrices, genotypes_detailed

    max_suppression = 1
    duration_of_protection = 0
    start_duration = 0
    max_Y = 0
    max_S = 0


    sex_breakpoint = 0
    for m in 1:length(genotypes_detailed)
        if genotypes_detailed[m][1] in X_chromosomes
            sex_breakpoint = m
            break
        end
    end

    # store_genotypes1 = deepcopy(genotype_vector1)
    # store_zygotes1 = zeros(length(genotype_vector1))
    # store_eggs1 = zeros(length(gametes_detailed))
    # store_sperm1 = zeros(length(gametes_detailed))

    # store_genotypes2 = deepcopy(genotype_vector2)
    # store_zygotes2 = zeros(length(genotype_vector2))
    # store_eggs2 = zeros(length(gametes_detailed))
    # store_sperm2 = zeros(length(gametes_detailed))

    Rm = current_Parameters["Rm"]
    theta = current_Parameters["theta"]

    #number of eggs ("f")
    f = (Rm * 2.0) / theta

    #alpha is a constant determining the density dependent mortality (half maximal)
    alpha = initial_size * f / (Rm - 1.0)

    for i in 1:t
        current_AllF = 0
        current_AllM = 0
        current_YLE_total = 0
        current_Shredder_total_all = 0


        post_migration1 = zeros(length(genotype_vector1))

        for j in 1:length(post_migration1)
            post_migration1[j] = ((genotype_vector1[j])
                                + (genotype_vector2[j] * migration_rate))
        end

        post_migration2 = zeros(length(genotype_vector2))

        for j in 1:length(post_migration2)
            post_migration2[j] = ((genotype_vector2[j] * (1 - migration_rate)))
        end

        genotype_vector1 = post_migration1
        genotype_vector2 = post_migration2

        #calculating prcesses in population 1
        genotype_vector1 = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector1, mutation_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector1, homing_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector1, editing_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector1, recombination_matrix_helper)

        sperm1 = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector1, sperm_matrix_helper)
        total_sperm1 = sum(sperm1)
        sperm1 = sperm1./total_sperm1

        eggs1 = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector1, egg_matrix_helper)
        eggs1 = eggs1 .* f

        # store_sperm1 = hcat(store_sperm1, sperm1)
        # store_eggs1 = hcat(store_eggs1, eggs1)

        genotype_vector1 = create_zygote_vector(sperm1, eggs1)

        # store_zygotes1 = hcat(store_zygotes1, genotype_vector1)

        #println("RUN ", i,)
        #println(sperm)
        #println(eggs)

        total_zygotes1 = sum(genotype_vector1)

        genotype_vector1 = genotype_vector1 .* (theta * alpha / (alpha + total_zygotes1))

        genotype_vector1 = current_matrices["selection_matrix"] .* genotype_vector1

        #calculating processes in population 2
        genotype_vector2 = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector2, mutation_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector2, homing_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector2, editing_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector2, recombination_matrix_helper)

        sperm2 = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector2, sperm_matrix_helper)
        total_sperm2 = sum(sperm2)
        sperm2 = sperm2./total_sperm2

        eggs2 = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector2, egg_matrix_helper)
        eggs2 = eggs2 .* f

        # store_sperm2 = hcat(store_sperm2, sperm2)
        # store_eggs2 = hcat(store_eggs2, eggs2)

        genotype_vector2 = create_zygote_vector(sperm2, eggs2)

        # store_zygotes2 = hcat(store_zygotes2, genotype_vector2)

        #println("RUN ", i,)
        #println(sperm)
        #println(eggs)

        total_zygotes2 = sum(genotype_vector2)

        genotype_vector2 = genotype_vector2 .* (theta * alpha / (alpha + total_zygotes2))

        genotype_vector2 = current_matrices["selection_matrix"] .* genotype_vector2


        # store_genotypes1 = hcat(store_genotypes1, genotype_vector1)
        #
        # store_genotypes2 = hcat(store_genotypes2, genotype_vector2)

        allfemales = 0
        for j in sex_breakpoint:length(genotypes_detailed)
            allfemales = allfemales + genotype_vector1[j]
        end

        for k in 1:length(genotypes_detailed)
            if (genotypes_detailed[k][1] in X_chromosomes)
                current_AllF += genotype_vector1[k]
            else
                current_AllM += genotype_vector1[k]
            end

            if (genotypes[k][1] == "Y2")
                current_YLE_total += genotype_vector1[k]
            end

            if ((genotypes_detailed[k][3] == "CD") && (genotypes_detailed[k][4] == "CD"))
                current_Shredder_total_all += (2 * genotype_vector1[k])
            elseif ((genotypes_detailed[k][3] == "CD") || (genotypes_detailed[k][4] == "CD"))
                current_Shredder_total_all += genotype_vector1[k]
            end
        end

        Y = current_YLE_total / current_AllM
        S = current_Shredder_total_all / (2*current_AllM + current_AllF)

        if Y > max_Y
            max_Y = Y
        end

        if S > max_S
            max_S = S
        end

        if allfemales < max_suppression
            max_suppression = allfemales
        end


        if (start_duration == 0) && (allfemales < threshold)
            start_duration = i
        elseif (start_duration != 0) && (allfemales > threshold)
            duration_of_protection = i - start_duration
            break
        end

    end

    # run1 = Dict(
    #     "genotypes" => store_genotypes1,
    #     # "zygotes" => store_zygotes1,
    #     # "eggs" => store_eggs1,
    #     # "sperm" => store_sperm1,
    # )
    #
    # run2 = Dict(
    #     "genotypes" => store_genotypes2,
    #     # "zygotes" => store_zygotes2,
    #     # "eggs" => store_eggs2,
    #     # "sperm" => store_sperm2,
    # )

    if (start_duration != 0) && (duration_of_protection == 0)
        duration_of_protection = t - start_duration
    end

    return max_suppression, duration_of_protection, max_Y, max_S

end

function duo_timecourse_parameterscan2_unidirectional_non_target(t, genotype_vector1, genotype_vector2;
                            rel_size = 1.0, migration_rate = 0.001,
                            initial_size = 1.0, threshold = 1e-2)
    global current_Parameters, current_matrices

    max_suppression = 1
    duration_of_protection = 0
    start_duration = 0
    max_Y = 0
    max_S = 0


    sex_breakpoint = 0
    for m in 1:length(genotypes_detailed)
        if genotypes_detailed[m][1] in X_chromosomes
            sex_breakpoint = m
            break
        end
    end

    # store_genotypes1 = deepcopy(genotype_vector1)
    # store_zygotes1 = zeros(length(genotype_vector1))
    # store_eggs1 = zeros(length(gametes_detailed))
    # store_sperm1 = zeros(length(gametes_detailed))

    # store_genotypes2 = deepcopy(genotype_vector2)
    # store_zygotes2 = zeros(length(genotype_vector2))
    # store_eggs2 = zeros(length(gametes_detailed))
    # store_sperm2 = zeros(length(gametes_detailed))

    Rm = current_Parameters["Rm"]
    theta = current_Parameters["theta"]

    #number of eggs ("f")
    f = (Rm * 2.0) / theta

    #alpha is a constant determining the density dependent mortality (half maximal)
    alpha = initial_size * f / (Rm - 1.0)

    for i in 1:t
        current_AllF = 0
        current_AllM = 0
        current_YLE_total = 0
        current_Shredder_total_all = 0

        post_migration1 = zeros(length(genotype_vector1))

        for j in 1:length(post_migration1)
            post_migration1[j] = ((genotype_vector1[j] * (1 - migration_rate)))

        end

        post_migration2 = zeros(length(genotype_vector2))

        for j in 1:length(post_migration2)
            post_migration2[j] = ((genotype_vector2[j])
                                + (genotype_vector1[j] * migration_rate))
        end

        genotype_vector1 = post_migration1
        genotype_vector2 = post_migration2

        #calculating prcesses in population 1
        genotype_vector1 = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector1, mutation_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector1, homing_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector1, editing_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector1, recombination_matrix_helper)

        sperm1 = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector1, sperm_matrix_helper)
        total_sperm1 = sum(sperm1)
        sperm1 = sperm1./total_sperm1

        eggs1 = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector1, egg_matrix_helper)
        eggs1 = eggs1 .* f

        # store_sperm1 = hcat(store_sperm1, sperm1)
        # store_eggs1 = hcat(store_eggs1, eggs1)

        genotype_vector1 = create_zygote_vector(sperm1, eggs1)

        # store_zygotes1 = hcat(store_zygotes1, genotype_vector1)

        #println("RUN ", i,)
        #println(sperm)
        #println(eggs)

        total_zygotes1 = sum(genotype_vector1)

        genotype_vector1 = genotype_vector1 .* (theta * alpha / (alpha + total_zygotes1))

        genotype_vector1 = current_matrices["selection_matrix"] .* genotype_vector1

        #calculating processes in population 2
        genotype_vector2 = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector2, mutation_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector2, homing_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector2, editing_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector2, recombination_matrix_helper)

        sperm2 = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector2, sperm_matrix_helper)
        total_sperm2 = sum(sperm2)
        sperm2 = sperm2./total_sperm2

        eggs2 = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector2, egg_matrix_helper)
        eggs2 = eggs2 .* f

        # store_sperm2 = hcat(store_sperm2, sperm2)
        # store_eggs2 = hcat(store_eggs2, eggs2)

        genotype_vector2 = create_zygote_vector(sperm2, eggs2)

        # store_zygotes2 = hcat(store_zygotes2, genotype_vector2)

        #println("RUN ", i,)
        #println(sperm)
        #println(eggs)

        total_zygotes2 = sum(genotype_vector2)

        genotype_vector2 = genotype_vector2 .* (theta * alpha / (alpha + total_zygotes2))

        genotype_vector2 = current_matrices["selection_matrix"] .* genotype_vector2


        # store_genotypes1 = hcat(store_genotypes1, genotype_vector1)
        #
        # store_genotypes2 = hcat(store_genotypes2, genotype_vector2)

        allfemales = 0
        for j in sex_breakpoint:length(genotypes_detailed)
            allfemales = allfemales + genotype_vector2[j]
        end

        if allfemales < max_suppression
            max_suppression = allfemales
        end

        for k in 1:length(genotypes_detailed)
            if (genotypes_detailed[k][1] in X_chromosomes)
                current_AllF += genotype_vector2[k]
            else
                current_AllM += genotype_vector2[k]
            end

            if (genotypes[k][1] == "Y2")
                current_YLE_total += genotype_vector2[k]
            end

            if ((genotypes_detailed[k][3] == "CD") && (genotypes_detailed[k][4] == "CD"))
                current_Shredder_total_all += (2 * genotype_vector2[k])
            elseif ((genotypes_detailed[k][3] == "CD") || (genotypes_detailed[k][4] == "CD"))
                current_Shredder_total_all += genotype_vector2[k]
            end
        end

        Y = current_YLE_total / current_AllM
        S = current_Shredder_total_all / (2*current_AllM + current_AllF)

        if Y > max_Y
            max_Y = Y
        end

        if S > max_S
            max_S = S
        end


        if (start_duration == 0) && (allfemales < threshold)
            start_duration = i
        elseif (start_duration != 0) && (allfemales > threshold)
            duration_of_protection = i - start_duration
            break
        end

    end

    # run1 = Dict(
    #     "genotypes" => store_genotypes1,
    #     # "zygotes" => store_zygotes1,
    #     # "eggs" => store_eggs1,
    #     # "sperm" => store_sperm1,
    # )
    #
    # run2 = Dict(
    #     "genotypes" => store_genotypes2,
    #     # "zygotes" => store_zygotes2,
    #     # "eggs" => store_eggs2,
    #     # "sperm" => store_sperm2,
    # )

    if (start_duration != 0) && (duration_of_protection == 0)
        duration_of_protection = t - start_duration
    end

    return max_suppression, duration_of_protection, max_Y, max_S

end


function duo_timecourse_parameterscan_extinction(t, genotype_vector1, genotype_vector2;
                            rel_size = 1.0, migration_rate = 0.001,
                            initial_size = 1.0, threshold = 1e-2, extinction_threshold = 1e-10)
    global current_Parameters, current_matrices

    max_suppression = 1
    duration_of_protection = 0
    start_duration = 0


    sex_breakpoint = 0
    for m in 1:length(genotypes_detailed)
        if genotypes_detailed[m][1] in X_chromosomes
            sex_breakpoint = m
            break
        end
    end

    # store_genotypes1 = deepcopy(genotype_vector1)
    # store_zygotes1 = zeros(length(genotype_vector1))
    # store_eggs1 = zeros(length(gametes_detailed))
    # store_sperm1 = zeros(length(gametes_detailed))

    # store_genotypes2 = deepcopy(genotype_vector2)
    # store_zygotes2 = zeros(length(genotype_vector2))
    # store_eggs2 = zeros(length(gametes_detailed))
    # store_sperm2 = zeros(length(gametes_detailed))

    Rm = current_Parameters["Rm"]
    theta = current_Parameters["theta"]

    #number of eggs ("f")
    f = (Rm * 2.0) / theta

    #alpha is a constant determining the density dependent mortality (half maximal)
    alpha = initial_size * f / (Rm - 1.0)

    for i in 1:t

        post_migration1 = zeros(length(genotype_vector1))

        for j in 1:length(post_migration1)
            post_migration1[j] = ((genotype_vector1[j] * (1 - migration_rate))
                                + (genotype_vector2[j] * migration_rate))
        end

        post_migration2 = zeros(length(genotype_vector2))

        for j in 1:length(post_migration2)
            post_migration2[j] = ((genotype_vector2[j] * (1 - migration_rate))
                                + (genotype_vector1[j] * migration_rate))
        end

        genotype_vector1 = post_migration1
        genotype_vector2 = post_migration2

        #calculating prcesses in population 1
        genotype_vector1 = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector1, mutation_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector1, homing_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector1, editing_matrix_helper)
        genotype_vector1 = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector1, recombination_matrix_helper)

        sperm1 = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector1, sperm_matrix_helper)
        total_sperm1 = sum(sperm1)
        sperm1 = sperm1./total_sperm1

        eggs1 = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector1, egg_matrix_helper)
        eggs1 = eggs1 .* f

        # store_sperm1 = hcat(store_sperm1, sperm1)
        # store_eggs1 = hcat(store_eggs1, eggs1)

        genotype_vector1 = create_zygote_vector(sperm1, eggs1)

        # store_zygotes1 = hcat(store_zygotes1, genotype_vector1)

        #println("RUN ", i,)
        #println(sperm)
        #println(eggs)

        total_zygotes1 = sum(genotype_vector1)

        genotype_vector1 = genotype_vector1 .* (theta * alpha / (alpha + total_zygotes1))

        genotype_vector1 = current_matrices["selection_matrix"] .* genotype_vector1

        #calculating processes in population 2
        genotype_vector2 = multiplication_optimised(current_matrices["mutation_matrix"], genotype_vector2, mutation_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["homing_matrix"], genotype_vector2, homing_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["editing_matrix"], genotype_vector2, editing_matrix_helper)
        genotype_vector2 = multiplication_optimised(current_matrices["recombination_matrix"], genotype_vector2, recombination_matrix_helper)

        sperm2 = multiplication_optimised(current_matrices["sperm_matrix"], genotype_vector2, sperm_matrix_helper)
        total_sperm2 = sum(sperm2)
        sperm2 = sperm2./total_sperm2

        eggs2 = multiplication_optimised(current_matrices["egg_matrix"], genotype_vector2, egg_matrix_helper)
        eggs2 = eggs2 .* f

        # store_sperm2 = hcat(store_sperm2, sperm2)
        # store_eggs2 = hcat(store_eggs2, eggs2)

        genotype_vector2 = create_zygote_vector(sperm2, eggs2)

        # store_zygotes2 = hcat(store_zygotes2, genotype_vector2)

        #println("RUN ", i,)
        #println(sperm)
        #println(eggs)

        total_zygotes2 = sum(genotype_vector2)

        genotype_vector2 = genotype_vector2 .* (theta * alpha / (alpha + total_zygotes2))

        genotype_vector2 = current_matrices["selection_matrix"] .* genotype_vector2


        # store_genotypes1 = hcat(store_genotypes1, genotype_vector1)
        #
        # store_genotypes2 = hcat(store_genotypes2, genotype_vector2)

        allfemales = 0
        for j in sex_breakpoint:length(genotypes_detailed)
            allfemales = allfemales + genotype_vector1[j]
        end

        if allfemales < max_suppression
            max_suppression = allfemales
        end


        if (start_duration == 0) && (allfemales < threshold)
            start_duration = i
        elseif (start_duration != 0) && (allfemales > threshold)
            duration_of_protection = i - start_duration
            break
        end

        if allfemales < extinction_threshold
            genotype_vector1 = zeros(length(genotype_vector1))
        end

    end

    # run1 = Dict(
    #     "genotypes" => store_genotypes1,
    #     # "zygotes" => store_zygotes1,
    #     # "eggs" => store_eggs1,
    #     # "sperm" => store_sperm1,
    # )
    #
    # run2 = Dict(
    #     "genotypes" => store_genotypes2,
    #     # "zygotes" => store_zygotes2,
    #     # "eggs" => store_eggs2,
    #     # "sperm" => store_sperm2,
    # )

    if (start_duration != 0) && (duration_of_protection == 0)
        duration_of_protection = t - start_duration
    end

    return max_suppression, duration_of_protection

end


#takes a process table and checks which terms in it are not zero
#it saves the non zero positions and returns them as a helper file
#which can then be passed to multiplication_optimised() to perform
# multiplications more efficiently
function create_helper(given_matrix)

    helper = []

    for i in 1:length(given_matrix[:,1])
        relevant_columns = []
        for j in 1:length(given_matrix[1,:])
            if given_matrix[i,j] != 0
                push!(relevant_columns, j)
            end
        end
        push!(helper, relevant_columns)
    end


    return helper
end


#takes a process table and converts all the values that are always
#unaffected by parameter values from type Symbol to type float.
function convert_matrix(given_matrix)


    remaining = []
    new_matrix = Array{Any}(undef,length(given_matrix[:,1]), length(given_matrix[1,:]))

    for i in 1:length(given_matrix[:,1])

        remaining_columns = []

        for j in 1:length(given_matrix[1,:])

            try
                new_matrix[i,j] = float(given_matrix[i,j])

            catch e
                push!(remaining_columns, j)
                new_matrix[i,j] = given_matrix[i,j]


            end

        end

        push!(remaining, remaining_columns)

    end

    return new_matrix, remaining

end


function one_dimensional_population_string(t, number_target, number_non_target,
        target_genotypes, non_target_genotypes; release_rate = 0.001,
        migration_rate = 0.001, initial_size = 1.0)

        Rm = current_Parameters["Rm"]
        theta = current_Parameters["theta"]

        #number of eggs ("f")
        f = (Rm * 2.0) / theta

        #alpha is a constant determining the density dependent mortality (half maximal)
        alpha = initial_size * f / (Rm - 1.0)

        #an array that contains a dictionary for every population in the string
        #only has 1 property: "genotypes"
        #so that the output can be used with other standard functions
        output = []

        #an array that also contains all populations and how they developed over time
        #but not yet in the form of a dictionary but as a simple array of tables
        output_genotypes = []

        #an array in which only the current makeup of every population is stored
        #itself also as an array of BigFloats
        current_populations = []


        #setting up a vector that contains the release population and
        #adding that as the first in the row
        release_pop = deepcopy(target_genotypes)
        release_pop[return_i(["AB" "ef" "CD" "CD"], genotypes_detailed)] = release_rate
        push!(current_populations, release_pop)

        #adding all the other target populations
        for i in 2:number_target
            push!(current_populations, deepcopy(target_genotypes))
        end

        #adding all the non-target populations
        for i in 1:number_non_target
            push!(current_populations, deepcopy(non_target_genotypes))
        end

        for i in 1:length(current_populations)
            push!(output_genotypes, current_populations[i])
        end


        for i in 1:t
            #first process every generation is that migration happens
            migrating = []

            #we calculate how many individuals are migrating in each population
            #we store the number of migrating individuals in "migrating"
            #and subtract it from the sending population

            for m in 1:length(current_populations)

                if ( (m == 1) || (m == length(current_populations)) )
                    push!(migrating, (current_populations[m].*migration_rate))
                    current_populations[m] = (current_populations[m].*(1 - migration_rate))
                else
                    push!(migrating, (current_populations[m].*(migration_rate)))
                    current_populations[m] = (current_populations[m].*(1 - (migration_rate *2)))
                end
            end

            #we add the immgrating individuals from left direction
            for m in 2:length(current_populations)
                for j in 1:length(current_populations[m])
                    current_populations[m][j] += migrating[(m - 1)][j]
                end
            end

            #we add the immgrating individuals from right direction
            for m in length(current_populations):-1:2
                for j in 1:length(current_populations[(m-1)])
                    current_populations[(m-1)][j] += migrating[m][j]
                end
            end

            #now we calculate the processes within each population
            for m in 1:length(current_populations)

                current_populations[m] = multiplication_optimised(current_matrices["mutation_matrix"], current_populations[m], mutation_matrix_helper)
                current_populations[m] = multiplication_optimised(current_matrices["homing_matrix"], current_populations[m], homing_matrix_helper)
                current_populations[m] = multiplication_optimised(current_matrices["editing_matrix"], current_populations[m], editing_matrix_helper)
                current_populations[m] = multiplication_optimised(current_matrices["recombination_matrix"], current_populations[m], recombination_matrix_helper)

                sperm = multiplication_optimised(current_matrices["sperm_matrix"], current_populations[m], sperm_matrix_helper)
                total_sperm = sum(sperm)
                sperm = sperm./total_sperm

                eggs = multiplication_optimised(current_matrices["egg_matrix"], current_populations[m], egg_matrix_helper)
                eggs = eggs .* f


                # println(sum(sperm))
                # println(sum(eggs))

                current_populations[m] = create_zygote_vector(sperm,eggs)

                total_zygotes = sum(current_populations[m])

                current_populations[m] = (current_populations[m] .* (theta * alpha / (alpha + total_zygotes)))


                current_populations[m] = (current_matrices["selection_matrix"] .* current_populations[m])

                output_genotypes[m] = hcat(output_genotypes[m], current_populations[m])

            end

        end

        for i in 1:length(current_populations)

            run = Dict(
                    "genotypes" => output_genotypes[i],
            )

            push!(output, run)

        end

        return output
end


function one_d_popstring_parameterscan(t, number_target, number_non_target,
        target_genotypes, non_target_genotypes; release_rate = 0.001,
        migration_rate = 0.001, initial_size = 1.0, threshold = 1e-2)

        Rm = current_Parameters["Rm"]
        theta = current_Parameters["theta"]

        #number of eggs ("f")
        f = (Rm * 2.0) / theta

        #alpha is a constant determining the density dependent mortality (half maximal)
        alpha = initial_size * f / (Rm - 1.0)

        # output = []
        # output_genotypes = []
        current_populations = []

        max_suppression = 1
        duration_of_protection = 0
        start_duration = 0


        #determine where the sex breakpoint is in the genotype list
        #for later calculations
        sex_breakpoint = 0
        for m in 1:length(genotypes_detailed)
            if genotypes_detailed[m][1] in X_chromosomes
                sex_breakpoint = m
                break
            end
        end


        #setting up a vector that contains the starting populations
        for i in 2:number_target
            push!(current_populations, target_genotypes)
        end

        target_genotypes[return_i(["AB" "ef" "CD" "CD"], genotypes_detailed)] = release_rate
        pushfirst!(current_populations, target_genotypes)

        for i in 1:number_non_target
            push!(current_populations, non_target_genotypes)
        end

        # for i in 1:length(current_populations)
        #     push!(output_genotypes, current_populations[i])
        #
        # end

        # println(typeof(current_populations))
        # println("_____")

        for i in 1:t

            #first migration happens
            migrating = []

            #we calculate how many individuals are lost in each population
            for m in 1:length(current_populations)

                if ( (m == 1) || (m == length(current_populations)) )
                    push!(migrating, (current_populations[m].*migration_rate))
                    current_populations[m] = (current_populations[m].*(1 - migration_rate))
                else
                    push!(migrating, (current_populations[m].*(migration_rate)))
                    current_populations[m] = (current_populations[m].*(1 - (migration_rate *2)))
                end
            end

            #we add the immgrating individuals from left direction
            for m in 2:length(current_populations)
                for j in 1:length(current_populations[m])
                    current_populations[m][j] += migrating[(m - 1)][j]
                end
            end

            #we add the immgrating individuals from right direction
            for m in length(current_populations):-1:2
                for j in 1:length(current_populations[(m-1)])
                    current_populations[(m-1)][j] += migrating[m][j]
                end
            end

            #now we calculate the processes within each population
            for m in 1:length(current_populations)

                current_populations[m] = multiplication_optimised(current_matrices["mutation_matrix"], current_populations[m], mutation_matrix_helper)
                current_populations[m] = multiplication_optimised(current_matrices["homing_matrix"], current_populations[m], homing_matrix_helper)
                current_populations[m] = multiplication_optimised(current_matrices["editing_matrix"], current_populations[m], editing_matrix_helper)
                current_populations[m] = multiplication_optimised(current_matrices["recombination_matrix"], current_populations[m], recombination_matrix_helper)

                sperm = multiplication_optimised(current_matrices["sperm_matrix"], current_populations[m], sperm_matrix_helper)
                total_sperm = sum(sperm)
                sperm = sperm./total_sperm

                eggs = multiplication_optimised(current_matrices["egg_matrix"], current_populations[m], egg_matrix_helper)
                eggs = eggs .* f


                # println(sum(sperm))
                # println(sum(eggs))

                current_populations[m] = create_zygote_vector(sperm,eggs)

                total_zygotes = sum(current_populations[m])

                current_populations[m] = (current_populations[m] .* (theta * alpha / (alpha + total_zygotes)))


                current_populations[m] = (current_matrices["selection_matrix"] .* current_populations[m])

                # output_genotypes[m] = hcat(output_genotypes[m], current_populations[m])

            end


            allfemales = 0

            #calculating the total number of females
            for m in 1:number_target
                allfemales += sum(current_populations[m][sex_breakpoint:end])
            end

            allfemales = allfemales/number_target


            if allfemales < max_suppression
                max_suppression = allfemales
            end


            if (start_duration == 0) && (allfemales < threshold)
                start_duration = i
            elseif (start_duration != 0) && (allfemales > threshold)
                duration_of_protection = i - start_duration
                break
            end

        end

        if (start_duration != 0) && (duration_of_protection == 0)
            duration_of_protection = t - start_duration
        end

        return max_suppression, duration_of_protection
end



function one_dimensional_population_string_with_extinction(t, number_target, number_non_target,
        target_genotypes, non_target_genotypes; release_rate = 0.001,
        migration_rate = 0.001, initial_size = 1.0, extinction_threshold = 1e-5)

        Rm = current_Parameters["Rm"]
        theta = current_Parameters["theta"]

        #number of eggs ("f")
        f = (Rm * 2.0) / theta

        #alpha is a constant determining the density dependent mortality (half maximal)
        alpha = initial_size * f / (Rm - 1.0)

        #determine where the sex breakpoint is in the genotype list
        #for later calculations
        sex_breakpoint = 0
        for m in 1:length(genotypes_detailed)
            if genotypes_detailed[m][1] in X_chromosomes
                sex_breakpoint = m
                break
            end
        end
        #an array that contains a dictionary for every population in the string
        #only has 1 property: "genotypes"
        #so that the output can be used with other standard functions
        output = []

        #an array that also contains all populations and how they developed over time
        #but not yet in the form of a dictionary but as a simple array of tables
        output_genotypes = []

        #an array in which only the current makeup of every population is stored
        #itself also as an array of BigFloats
        current_populations = []


        #setting up a vector that contains the release population and
        #adding that as the first in the row
        release_pop = deepcopy(target_genotypes)
        release_pop[return_i(["AB" "ef" "CD" "CD"], genotypes_detailed)] = release_rate
        push!(current_populations, release_pop)

        #adding all the other target populations
        for i in 2:number_target
            push!(current_populations, deepcopy(target_genotypes))
        end

        #adding all the non-target populations
        for i in 1:number_non_target
            push!(current_populations, deepcopy(non_target_genotypes))
        end

        for i in 1:length(current_populations)
            push!(output_genotypes, current_populations[i])
        end


        for i in 1:t

            #first migration happens
            migrating = []

            #we calculate how many individuals are lost in each population
            for m in 1:length(current_populations)

                if ( (m == 1) || (m == length(current_populations)) )
                    push!(migrating, (current_populations[m].*migration_rate))
                    current_populations[m] = (current_populations[m].*(1 - migration_rate))
                else
                    push!(migrating, (current_populations[m].*(migration_rate)))
                    current_populations[m] = (current_populations[m].*(1 - (migration_rate *2)))
                end
            end

            #we add the immgrating individuals from left direction
            for m in 2:length(current_populations)
                for j in 1:length(current_populations[m])
                    current_populations[m][j] += migrating[(m - 1)][j]
                end
            end

            #we add the immgrating individuals from right direction
            for m in length(current_populations):-1:2
                for j in 1:length(current_populations[(m-1)])
                    current_populations[(m-1)][j] += migrating[m][j]
                end
            end

            #now we calculate the processes within each population
            for m in 1:length(current_populations)

                current_populations[m] = multiplication_optimised(current_matrices["mutation_matrix"], current_populations[m], mutation_matrix_helper)
                current_populations[m] = multiplication_optimised(current_matrices["homing_matrix"], current_populations[m], homing_matrix_helper)
                current_populations[m] = multiplication_optimised(current_matrices["editing_matrix"], current_populations[m], editing_matrix_helper)
                current_populations[m] = multiplication_optimised(current_matrices["recombination_matrix"], current_populations[m], recombination_matrix_helper)

                sperm = multiplication_optimised(current_matrices["sperm_matrix"], current_populations[m], sperm_matrix_helper)
                total_sperm = sum(sperm)
                sperm = sperm./total_sperm

                eggs = multiplication_optimised(current_matrices["egg_matrix"], current_populations[m], egg_matrix_helper)
                eggs = eggs .* f


                # println(sum(sperm))
                # println(sum(eggs))

                current_populations[m] = create_zygote_vector(sperm,eggs)

                total_zygotes = sum(current_populations[m])

                current_populations[m] = (current_populations[m] .* (theta * alpha / (alpha + total_zygotes)))


                current_populations[m] = (current_matrices["selection_matrix"] .* current_populations[m])

                all_females = sum(current_populations[m][sex_breakpoint:end])

                if all_females < extinction_threshold
                    current_populations[m] = zeros(length(current_populations[m]))
                end

                output_genotypes[m] = hcat(output_genotypes[m], current_populations[m])

            end

        end

        for i in 1:length(current_populations)

            run = Dict(
                    "genotypes" => output_genotypes[i],
            )

            push!(output, run)

        end

        return output
end



function one_d_popstring_extinction_parameterscan(t, number_target, number_non_target,
        target_genotypes, non_target_genotypes; release_rate = 0.001,
        migration_rate = 0.001, initial_size = 1.0, threshold = 1e-2, extinction_threshold = 1e-10)

        Rm = current_Parameters["Rm"]
        theta = current_Parameters["theta"]

        #number of eggs ("f")
        f = (Rm * 2.0) / theta

        #alpha is a constant determining the density dependent mortality (half maximal)
        alpha = initial_size * f / (Rm - 1.0)


        max_suppression = 1
        duration_of_protection = 0
        start_duration = 0


        #determine where the sex breakpoint is in the genotype list
        #for later calculations
        sex_breakpoint = 0
        for m in 1:length(genotypes_detailed)
            if genotypes_detailed[m][1] in X_chromosomes
                sex_breakpoint = m
                break
            end
        end


        #an array in which only the current makeup of every population is stored
        #itself also as an array of BigFloats
        current_populations = []


        #setting up a vector that contains the release population and
        #adding that as the first in the row
        release_pop = deepcopy(target_genotypes)
        release_pop[return_i(["AB" "ef" "CD" "CD"], genotypes_detailed)] = release_rate
        push!(current_populations, release_pop)

        #adding all the other target populations
        for i in 2:number_target
            push!(current_populations, deepcopy(target_genotypes))
        end

        #adding all the non-target populations
        for i in 1:number_non_target
            push!(current_populations, deepcopy(non_target_genotypes))
        end

        # for i in 1:length(current_populations)
        #     push!(output_genotypes, current_populations[i])
        #
        # end

        # println(typeof(current_populations))
        # println("_____")

        for i in 1:t

            #first migration happens
            migrating = []

            #we calculate how many individuals are lost in each population
            for m in 1:length(current_populations)

                if ( (m == 1) || (m == length(current_populations)) )
                    push!(migrating, (current_populations[m].*migration_rate))
                    current_populations[m] = (current_populations[m].*(1 - migration_rate))
                else
                    push!(migrating, (current_populations[m].*(migration_rate)))
                    current_populations[m] = (current_populations[m].*(1 - (migration_rate *2)))
                end
            end

            #we add the immgrating individuals from left direction
            for m in 2:length(current_populations)
                for j in 1:length(current_populations[m])
                    current_populations[m][j] += migrating[(m - 1)][j]
                end
            end

            #we add the immgrating individuals from right direction
            for m in length(current_populations):-1:2
                for j in 1:length(current_populations[(m-1)])
                    current_populations[(m-1)][j] += migrating[m][j]
                end
            end

            #now we calculate the processes within each population
            for m in 1:length(current_populations)

                current_populations[m] = multiplication_optimised(current_matrices["mutation_matrix"], current_populations[m], mutation_matrix_helper)
                current_populations[m] = multiplication_optimised(current_matrices["homing_matrix"], current_populations[m], homing_matrix_helper)
                current_populations[m] = multiplication_optimised(current_matrices["editing_matrix"], current_populations[m], editing_matrix_helper)
                current_populations[m] = multiplication_optimised(current_matrices["recombination_matrix"], current_populations[m], recombination_matrix_helper)

                sperm = multiplication_optimised(current_matrices["sperm_matrix"], current_populations[m], sperm_matrix_helper)
                total_sperm = sum(sperm)
                sperm = sperm./total_sperm

                eggs = multiplication_optimised(current_matrices["egg_matrix"], current_populations[m], egg_matrix_helper)
                eggs = eggs .* f


                # println(sum(sperm))
                # println(sum(eggs))

                current_populations[m] = create_zygote_vector(sperm,eggs)

                total_zygotes = sum(current_populations[m])

                current_populations[m] = (current_populations[m] .* (theta * alpha / (alpha + total_zygotes)))


                current_populations[m] = (current_matrices["selection_matrix"] .* current_populations[m])

                allfemales = sum(current_populations[m][sex_breakpoint:end])

                if allfemales < extinction_threshold
                    current_populations[m] = zeros(length(current_populations[m]))
                end

                # output_genotypes[m] = hcat(output_genotypes[m], current_populations[m])

            end


            allfemales = 0

            #calculating the total number of females
            for m in 1:number_target
                allfemales += sum(current_populations[m][sex_breakpoint:end])
            end

            allfemales = allfemales/number_target


            if allfemales < max_suppression
                max_suppression = allfemales
            end


            if (start_duration == 0) && (allfemales < threshold)
                start_duration = i
            elseif (start_duration != 0) && (allfemales > threshold)
                duration_of_protection = i - start_duration
                break
            end

        end

        if (start_duration != 0) && (duration_of_protection == 0)
            duration_of_protection = t - start_duration
        end

        return max_suppression, duration_of_protection
end


function sensitivity_plot(scan; ymin = 1e-10, kwargs...)

    # plot(scan[2:end,1], scan[2:end,2], lw = 4, label = "min pop size",
    #     dpi = 300, size = (500,300))
    # plot!(scan[2:end,1], scan[2:end,4], lw = 4, label = "transgenic Y",)
    # plot!(scan[2:end,1], scan[2:end,5], lw = 4, label = "transgenic autosomes",)
    # plot!(legend= false)
    #
    # plot!(twinx(), scan[2:end,1], scan[2:end,3], lw = 4, label = "duration below 1%",
    #     dpi = 300, size = (500,300), colour = "brown", legend = false)

    plot(scan[2:end,1], (scan[2:end,2]), lw = 4, label = "min pop size", colour=:black,
        linestyle=:dot, dpi = 300, size = (325,150), legend = false, ylims = (ymin, 5), yaxis =:log)
    plot!(ytickfontsize = 10, xtickfontsize = 10)
    # plot!(twinx(), scan[2:end,1], scan[2:end,3], lw = 4, label = "duration below 1%",
    #     dpi = 300, size = (325,150), colour =:black, legend = false, ylims = (-2, 152))
    plot!(grid = false)
    # plot!(ytickfontsize = 10, xtickfontsize = 10)
end


function sensitivity_plot2(scan1, scan2; given_size = (325,200), ymin = 1e-10, kwargs...)


    plot(scan1[2:end,1], (scan1[2:end,2]), lw = 4, label = "min pop size", colour=:black,
        linestyle=:dot,
        dpi = 300, size = given_size, legend = false, ylims = (ymin, 5), yaxis =:log)
    plot!(ytickfontsize = 10, xtickfontsize = 10, grid = false)
    plot!(scan1[2:end,1], (scan2[2:end,2]), lw = 4, colour=:black, dpi = 300, legend = false, ylims = (ymin, 5), yaxis =:log, linestyle=:dash)
    #plot!(scan1[2:end,1], (scan3[2:end,2]), lw = 4,  colour=:blue, dpi = 300, size = (325,150), legend = false, ylims = (ymin, 5), yaxis =:log)

end


function sensitivity_plot3(scan1, scan2, scan3, scan4, scan5, scan6;
            given_size = (300,180), ymin = 1e-10, kwargs...)


    plot(scan1[2:end,1], (scan1[2:end,2]), lw = 4, label = "min pop size", colour= colorant"#000000",
        linestyle=:dot,
        dpi = 300, size = given_size, legend = false, ylims = (ymin, 5), yaxis =:log)
    plot!(ytickfontsize = 10, xtickfontsize = 10, grid = false)
    plot!(scan1[2:end,1], (scan2[2:end,2]), lw = 4, colour= colorant"#000000", dpi = 300, legend = false, ylims = (ymin, 5), yaxis =:log, linestyle=:dash)
    plot!(scan1[2:end,1], (scan3[2:end,2]), lw = 4, colour= colorant"#F0DC00", dpi = 300, legend = false, ylims = (ymin, 5), yaxis =:log, linestyle=:dot)
    plot!(scan1[2:end,1], (scan4[2:end,2]), lw = 4, colour= colorant"#F0DC00", dpi = 300, legend = false, ylims = (ymin, 5), yaxis =:log, linestyle=:dash)
    plot!(scan1[2:end,1], (scan5[2:end,2]), lw = 4, colour= colorant"#FF7200", dpi = 300, legend = false, ylims = (ymin, 5), yaxis =:log, linestyle=:dot)
    plot!(scan1[2:end,1], (scan6[2:end,2]), lw = 4, colour= colorant"#FF7200", dpi = 300, legend = false, ylims = (ymin, 5), yaxis =:log, linestyle=:dash)
    plot!(ytickfontsize = 12, xtickfontsize = 12, dpi = 300, size = given_size)

end


function return_stats(given_run; lower_sup = 1e-10, upper_sup = 0.95, max_transgenic = 0.1)

    lower = 0.0
    time_limit = 0.0
    upper = 100.0

    sup_column = findfirst(x->x=="max_sup", given_run[1,:])
    Y_column = findfirst(x->x=="Y", given_run[1,:])
    S_column = findfirst(x->x=="S", given_run[1,:])



    for i in 2:length(given_run[:,1])
        if (given_run[i,sup_column] < lower_sup)
            lower = i
        end

        if (given_run[i,4] > 50) && (time_limit == 0.0)
            time_limit = i
        end

        if (given_run[i, sup_column]> upper_sup) && (given_run[i,Y_column]< max_transgenic) && (given_run[i,S_column]< max_transgenic)
            upper = i
            break
        end

    end

    println("lower bounday: ", given_run[lower, 1])
    println("upper boundary: ", given_run[upper,1])
    println("time to reach 1e-2 within 50 generations is at ", given_run[time_limit,1])
end


function return_min_population_size_and_max_transgenic_freq(input)
    global current_matrices
    output_min_population_size = BigFloat("1")
    output_max_transgenic_Y = BigFloat("0")
    output_max_transgenic_S = BigFloat("0")
    output_min_female_proportion = BigFloat("0.5")
    output_max_load = BigFloat("0")

    for t in 1:length(input["genotypes"][1,:])

        allfemales = BigFloat("0")
        allmales = BigFloat("0")
        cur_transgenic_Y = BigFloat("0")
        cur_transgenic_S = BigFloat("0")
        cur_female_proportion = BigFloat("0.5")
        cur_load = BigFloat("0")
        current_AllF_zygotes = BigFloat("0")
        current_AllF_x_selection_zygotes = BigFloat("0")

        genotype_vector = input["genotypes"][:,t]

        for k in 1:length(genotype_vector)

            if genotypes_detailed[k][1] in X_chromosomes
                allfemales += genotype_vector[k]
                current_AllF_zygotes += input["zygotes"][k,t]
                current_AllF_x_selection_zygotes += input["zygotes"][k,t] * current_matrices["selection_matrix"][k]
            else
                allmales += genotype_vector[k]
            end

            if occursin("AB", genotypes_detailed[k][1])
                cur_transgenic_Y +=  genotype_vector[k]
            elseif occursin("A", genotypes_detailed[k][1])
                cur_transgenic_Y +=  genotype_vector[k]
            elseif occursin("B", genotypes_detailed[k][1])
                cur_transgenic_Y +=  genotype_vector[k]
            end


            if ((occursin("C", genotypes_detailed[k][3]) || occursin("D", genotypes_detailed[k][3]))
                &&  (occursin("C", genotypes_detailed[k][4])) || occursin("D", genotypes_detailed[k][4]))
                cur_transgenic_S += 2 * genotype_vector[k]

            elseif ((occursin("C", genotypes_detailed[k][3]) || occursin("D", genotypes_detailed[k][3]))
                ||  (occursin("C", genotypes_detailed[k][4])) || occursin("D", genotypes_detailed[k][4]))
                cur_transgenic_S +=  genotype_vector[k]
            end

        end

        cur_transgenic_Y = cur_transgenic_Y / allmales
        cur_transgenic_S = cur_transgenic_S / (allmales * 2 + allfemales * 2)
        cur_female_proportion = allfemales/(allmales + allfemales)

        cur_load = 1 - 2 * (current_AllF_zygotes/sum(
                    input["zygotes"][:,t])) * ((current_AllF_x_selection_zygotes) / current_AllF_zygotes)

        if allfemales < output_min_population_size
            output_min_population_size = allfemales
        end

        if cur_transgenic_Y > output_max_transgenic_Y
            output_max_transgenic_Y = cur_transgenic_Y
        end

        if cur_transgenic_S > output_max_transgenic_S
            output_max_transgenic_S = cur_transgenic_S
        end

        if cur_female_proportion < output_min_female_proportion
            output_min_female_proportion = cur_female_proportion
        end

        if cur_load > output_max_load
            output_max_load = cur_load
        end

    end



    println("Minimum population size is ", output_min_population_size)
    println("Max transgenic Y is ", output_max_transgenic_Y)
    println("Max transgenic Autosome is ", output_max_transgenic_S)
    println("Min female proportion is ", output_min_female_proportion)
    println("Max genetic load is ", output_max_load)

end



function return_geographical_spread(input)
    #takes a run file of a one dimensional population string
    #determines for each population what the min population size and max frequency
    #of transgenic elements was
    #and returns a table of those values so that they can be plotted or processed

    output = Any["number" "sup" "Y" "S"]

    #iterate through all populations
    for i in 1:length(input)

        max_suppression = BigFloat("1")
        transgenic_Y = BigFloat("0")
        transgenic_S = BigFloat("0")

        #iterate through all generations in each population
        for k in 1:length(input[i]["genotypes"][1,:])

            genotype_vector = input[i]["genotypes"][:,k]

            allfemales = BigFloat("0")
            allmales = BigFloat("0")
            cur_transgenic_Y = BigFloat("0")
            cur_transgenic_S = BigFloat("0")

            #iterate through all individuals in each generation
            for j in 1:length(genotype_vector)


                if genotypes_detailed[j][1] in X_chromosomes
                    allfemales += genotype_vector[j]
                else
                    allmales += genotype_vector[j]
                end

                if occursin("AB", genotypes_detailed[j][1])
                    cur_transgenic_Y +=  genotype_vector[j]
                elseif occursin("A", genotypes_detailed[j][1])
                    cur_transgenic_Y +=  genotype_vector[j]
                elseif occursin("B", genotypes_detailed[j][1])
                    cur_transgenic_Y +=  genotype_vector[j]
                end


                if ((occursin("C", genotypes_detailed[j][3]) || occursin("D", genotypes_detailed[j][3]))
                    &&  (occursin("C", genotypes_detailed[j][4])) || occursin("D", genotypes_detailed[j][4]))
                    cur_transgenic_S += 2 * genotype_vector[j]

                elseif ((occursin("C", genotypes_detailed[j][3]) || occursin("D", genotypes_detailed[j][3]))
                    ||  (occursin("C", genotypes_detailed[j][4])) || occursin("D", genotypes_detailed[j][4]))
                    cur_transgenic_S +=  genotype_vector[j]
                end

            end

            cur_transgenic_Y = cur_transgenic_Y / allmales
            cur_transgenic_S = cur_transgenic_S / (allmales * 2 + allfemales * 2)


            if allfemales < max_suppression
                max_suppression = allfemales
            end

            if cur_transgenic_Y > transgenic_Y
                transgenic_Y = cur_transgenic_Y
            end

            if cur_transgenic_S > transgenic_S
                transgenic_S = cur_transgenic_S
            end

        end

        add = [i max_suppression transgenic_Y transgenic_S]
        output = vcat(output, add)

    end

    return output

end

# ______________________________________________________________________________
# initialise section
# ______________________________________________________________________________

current_Parameters = Dict();
current_matrices = Dict();

selection_matrix = create_selection_matrix(genotypes_detailed)
sperm_matrix, egg_matrix = create_gamete_matrix(genotypes_detailed, gametes_detailed)
editing_matrix = create_editing_matrix(genotypes_detailed)
homing_matrix = create_homing_matrix(genotypes_detailed)
mutation_matrix = create_mutation_matrix(genotypes_detailed)
recombination_matrix = create_recombination_matrix(genotypes_detailed)


sperm_matrix_helper = create_helper(sperm_matrix)
egg_matrix_helper = create_helper(egg_matrix)
editing_matrix_helper = create_helper(editing_matrix)
homing_matrix_helper = create_helper(homing_matrix)
mutation_matrix_helper = create_helper(mutation_matrix)
recombination_matrix_helper = create_helper(recombination_matrix)

sperm_matrix_num, sperm_symbols = convert_matrix(sperm_matrix)
egg_matrix_num, egg_symbols = convert_matrix(egg_matrix)
editing_matrix_num, editing_symbols = convert_matrix(editing_matrix)
homing_matrix_num, homing_symbols = convert_matrix(homing_matrix)
mutation_matrix_num, mutation_symbols = convert_matrix(mutation_matrix)
recombination_matrix_num, recombination_symbols = convert_matrix(recombination_matrix)


#here we define the ideal baseline parameter set
Parameters_set_standard = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1.0"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("1"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("1"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0"),

    "s_a" => BigFloat("0"),
    "s_b" => BigFloat("0"),
    "s_c" => BigFloat("0"),
    "s_d" => BigFloat("0"),
    "s_e" => BigFloat("0"),
    "h_e" => BigFloat("1"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.953"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("0"),

    #background mutation rate of all components
    "m_2" => BigFloat("0"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.947368"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("1"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("6"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("0"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("0"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)


#here we define the realistic baseline parameter set
Parameters_set_realistic = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1.0"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("0.99"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("0.95"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0.01"),

    "s_a" => BigFloat("0.01"),
    "s_b" => BigFloat("0.005"),
    "s_c" => BigFloat("0.01"),
    "s_d" => BigFloat("0.01"),
    "s_e" => BigFloat("0.01"),
    "h_e" => BigFloat("0.99"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.953"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("1e-3"),

    #background mutation rate of all components
    "m_2" => BigFloat("1e-6"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.947368"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("0.99"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("6"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("1e-6"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("1e-3"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)

#here we define a mixed ideal / realistic baseline parameter set
Parameters_set_realistic2 = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1.0"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("0.99"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("0.95"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0.0"),

    "s_a" => BigFloat("0.0"),
    "s_b" => BigFloat("0.00"),
    "s_c" => BigFloat("0.0"),
    "s_d" => BigFloat("0.0"),
    "s_e" => BigFloat("0.0"),
    "h_e" => BigFloat("1"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.953"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("1e-3"),

    #background mutation rate of all components
    "m_2" => BigFloat("1e-6"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.947368"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("0.99"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("6"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("1e-6"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("1e-3"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)



#here we define the realistic baseline parameter set
Parameters_set_realistic = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1.0"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("0.99"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("0.95"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0.01"),

    "s_a" => BigFloat("0.01"),
    "s_b" => BigFloat("0.005"),
    "s_c" => BigFloat("0.01"),
    "s_d" => BigFloat("0.01"),
    "s_e" => BigFloat("0.01"),
    "h_e" => BigFloat("0.99"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.953"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("1e-3"),

    #background mutation rate of all components
    "m_2" => BigFloat("1e-6"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.947368"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("0.99"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("6"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("1e-6"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("1e-3"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)


#another parameter set
Parameters_set_sensitivity1 = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("1"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("1"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0"),

    "s_a" => BigFloat("0"),
    "s_b" => BigFloat("0"),
    "s_c" => BigFloat("0"),
    "s_d" => BigFloat("0.01"),
    "s_e" => BigFloat("0.01"),
    "h_e" => BigFloat("1"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.95"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("1e-3"),

    #background mutation rate of all components
    "m_2" => BigFloat("1e-6"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.90"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("1.0"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("6"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("0"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("0"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)

#another parameter set
Parameters_set_sensitivity2 = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("1"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("1"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0"),

    "s_a" => BigFloat("0"),
    "s_b" => BigFloat("0"),
    "s_c" => BigFloat("0"),
    "s_d" => BigFloat("0.01"),
    "s_e" => BigFloat("0.01"),
    "h_e" => BigFloat("1"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.95"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("1e-3"),

    #background mutation rate of all components
    "m_2" => BigFloat("1e-6"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.90"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("1.0"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("6"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("0"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("0.05"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)

#another parameter set
Parameters_set_sensitivity3 = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("1"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("1"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0"),

    "s_a" => BigFloat("0"),
    "s_b" => BigFloat("0"),
    "s_c" => BigFloat("0"),
    "s_d" => BigFloat("0.01"),
    "s_e" => BigFloat("0.01"),
    "h_e" => BigFloat("1"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.95"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("1e-3"),

    #background mutation rate of all components
    "m_2" => BigFloat("1e-6"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.90"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("1.0"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("6"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("0"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("0.6"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)


#another parameter set
Parameters_set_sensitivity4_lowRM = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("1"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("1"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0"),

    "s_a" => BigFloat("0"),
    "s_b" => BigFloat("0"),
    "s_c" => BigFloat("0"),
    "s_d" => BigFloat("0.01"),
    "s_e" => BigFloat("0.01"),
    "h_e" => BigFloat("1"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.95"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("1e-3"),

    #background mutation rate of all components
    "m_2" => BigFloat("1e-6"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.90"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("1.0"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("4"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("0"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("0.05"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)

#another parameter set
Parameters_set_sensitivity5_lowRM = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("1"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("1"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0"),

    "s_a" => BigFloat("0"),
    "s_b" => BigFloat("0"),
    "s_c" => BigFloat("0"),
    "s_d" => BigFloat("0.01"),
    "s_e" => BigFloat("0.01"),
    "h_e" => BigFloat("1"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.95"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("1e-3"),

    #background mutation rate of all components
    "m_2" => BigFloat("1e-6"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.90"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("1.0"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("4"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("0"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("0.6"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)


#another parameter set
Parameters_set_sensitivity6_highRM = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("1"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("1"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0"),

    "s_a" => BigFloat("0"),
    "s_b" => BigFloat("0"),
    "s_c" => BigFloat("0"),
    "s_d" => BigFloat("0.01"),
    "s_e" => BigFloat("0.01"),
    "h_e" => BigFloat("1"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.95"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("1e-3"),

    #background mutation rate of all components
    "m_2" => BigFloat("1e-6"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.90"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("1.0"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("8"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("0"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("0.05"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)

#another parameter set
Parameters_set_sensitivity7_highRM = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("1"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("1"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0"),

    "s_a" => BigFloat("0"),
    "s_b" => BigFloat("0"),
    "s_c" => BigFloat("0"),
    "s_d" => BigFloat("0.01"),
    "s_e" => BigFloat("0.01"),
    "h_e" => BigFloat("1"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0.95"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0.95"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("1e-3"),

    #background mutation rate of all components
    "m_2" => BigFloat("1e-6"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0.90"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("1.0"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("8"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("0"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("0.6"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)


#wildtye parameters to calculate genetic load
Parameters_set_wildtype = Dict(

    #conditionality factor, describing whether X-shredding only happens in presence
    #of YLE (=0) or in presence of all Y chromosomes (=1)
    "c" => BigFloat("1"),

    # selection against carrying edited X gene (homozygous, females)
    "s_f" => BigFloat("0"),

    # dominance coefficient for carrying edited X gene (heterozygous, females)
    "h_f" => BigFloat("1"),

    # selection against carrying edited X gene (males)
    "s_m" => BigFloat("0"),

    "s_a" => BigFloat("0"),
    "s_b" => BigFloat("0"),
    "s_c" => BigFloat("0"),
    "s_d" => BigFloat("0"),
    "s_e" => BigFloat("0"),
    "h_e" => BigFloat("1"),

    #efficiency of editing (e_e)
    "e_e" => BigFloat("0"),

    #efficiency of homing (e_c)
    "e_h" => BigFloat("0"),

    #mutation rate during homing (m)
    "m_1" => BigFloat("0"),

    #background mutation rate of all components
    "m_2" => BigFloat("0"),

    #efficiency of X-shredding
    "e_s" => BigFloat("0"),

    #dominance coefficient of X-shredding activity
    "h_e2" => BigFloat("1.0"),

    #intrinsic rate of population increase
    "Rm" => BigFloat("6"),

    #density independent survival rate
    "theta" => BigFloat("0.1"),

    #rate with which editing resistance occurs during editing process
    "er_1" => BigFloat("0"),

    #rate with which shredding resistance occurs during shredding process
    "er_2" => BigFloat("0"),

    #rate with which homing resistance occurs during the homing process
    "er_3" => BigFloat("0"),

    #recombination rate (X-linked loci;Editing gene and shredding target site)
    "r" => BigFloat("0.5"),
)

#we calculate tables for wild-type conditions (no editing, no homing, no fitness effects...)
apply_parameters_set(Parameters_set_wildtype)
wildtype_matrices = deepcopy(current_matrices)

#we apply the baseline parameter values to our tables
#so that we can run timecourses
apply_parameters_set(Parameters_set_standard)


#we perform a plausibility check on the symbolic terms in each row of our
#process tables
plausibility_check(homing_matrix, mutation_matrix, egg_matrix, sperm_matrix,
                    editing_matrix, recombination_matrix)

# plausibility_check2(current_matrices["homing_matrix"], current_matrices["mutation_matrix"],
#                     current_matrices["egg_matrix"], current_matrices["sperm_matrix"],
#                     current_matrices["editing_matrix"], current_matrices["recombination_matrix"])



genotype_standard = generate_population();
genotype_standard[return_i(["AB" "ef" "CD" "cd"], genotypes_detailed)] = BigFloat("0.001");

#non_target = generate_population(r2 = 0.5, r3 = 0.5);
