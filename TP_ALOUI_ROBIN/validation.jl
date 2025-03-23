using DelimitedFiles
using LinearAlgebra
using SparseArrays
using Plots

function f1(x) # f qui a en entree x un nombre reel
    return 1.0 # sortie f(x) pre-remplit une valeur arbitraire
end

function f2(x) # f qui a en entree x un nombre reel
    return sin(x) # sortie f(x) pre-remplit une valeur arbitraire
end

function maillage(taille,NE,nom_fichier)
    out=open(nom_fichier,"w")
    NN=NE+1
    # Nombre de noeuds NN et Nombre d'elements NE
    println(out,NN," ",NE)
    # Numero des noeud (tabulation) et coordonnees des noeud
    h=taille/NE
    for i in 0:(NN-1) # on va jusqu'a NN-1 car on fait i+1 dans le print
        println(out,i+1," ",i*h)  # Indices des nœuds commencent à 1
    end
    # Numeros des elements et liste des noeuds de l'element
    for j in 1:NE
        println(out,j," ",j," ",j+1)
    end
    # Nombre de noeuds ND avec des conditions limites de Dirichlet
    println(out,2)
    # Numeros des noeuds avec conditions de Dirichlet
    println(out,1)
    println(out,NN)
    close(out)
end

function lecture_maillage(nom_fichier)
    #les positions des noeuds, les elements et les noeuds sur conditions de Dirichlet
    data=readdlm(nom_fichier) # le résultat est une matrice contenant toutes les valeurs du fichier
    # le nombre de noeuds
    NN=Int(data[1,1])
    # le nombre d'elements
    NE=Int(data[1,2])
    # les positions des noeuds
    liste_noeuds=data[2:NN+1,2]
    # liste des noeuds de l'element
    liste_elements=data[NN+2:NN+1+NE,2:3]
    # les noeuds sur conditions de Dirichlet
    noeud_dirichlet=data[NN+1+NE+2:NN+NE+4,1]
    return liste_noeuds,liste_elements,noeud_dirichlet
end

function tables_elementaires(num_element,liste_noeuds,liste_elements,f)
    # On extrait les indices des nœuds de l'élément
    noeud_1=Int(liste_elements[num_element,1])
    noeud_2=Int(liste_elements[num_element,2])
    # On recupere les coordonnées des deux noeuds
    coord_1=liste_noeuds[Int(noeud_1)]
    coord_2=liste_noeuds[Int(noeud_2)]
    # On calcul la longueur de l'élément
    h=coord_2-coord_1
     # On définit de la matrice de rigidité
    Ke=(1.0/h)*[1 -1;-1 1]
    # On évalue f(x) aux extrémités de l'élément
    f_coord1=f(coord_1)
    f_coord2=f(coord_2)
    # On calcul le second membre avec la formule du Trapèze
    Fe=(h/2)*[f_coord1;f_coord2]
    return Ke,Fe
end


function assembler_matrices(liste_noeuds,liste_elements,f)
    # On récupère le nombre de ligne qui correspodn au nombre d'element
    NE=size(liste_elements,1)  # Probleme quand on utilise length on a le nombre d'elment total
    # On récupère le nombre total de noeuds
    NN=length(liste_noeuds)
    # On initialise la matrice de rigidité globale
    A =zeros(NN, NN)
    # On initialise le vecteur du second membre global
    b=zeros(NN)
    # On assemblage les matrices élémentaires
    for i in 1:NE
        # On récupère les deux noeuds de l'elment actuel
        noeud1,noeud2=liste_elements[i,:]
        # On calcul les tables elementaires
        Ke,Fe=tables_elementaires(i,liste_noeuds,liste_elements,f)
        # On assemblage la matrice de rigidité globale
        A[noeud1,noeud1]+=Ke[1, 1]
        A[noeud2,noeud1]+=Ke[2, 1]
        A[noeud1,noeud2]+=Ke[1, 2]
        A[noeud2,noeud2]+=Ke[2, 2]

        # On fait l'assemblage dans le vecteur second membre
        b[noeud1]+=Fe[1]
        b[noeud2]+=Fe[2]
    end
    return A,b
end

function annuler_ligne_colonne(mat,indices)
    # On assigne chaque élément individuellement sans créer de nouvelle copie
    mat[indices,:].=0.0  # On met à zéro la ligne idx
    mat[:,indices].=0.0  # On met à zéro la colonne idx
end


function imposer_cl_dirichlet(A_in,b_in,noeud_dirichlet)
    # On fait une copie des données d'entrée pour préserver les originaux
    A_out=deepcopy(A_in)
    b_out=deepcopy(b_in)
    # On récupère la taille du système en récupérant le nbr de ligne (1)
    n=size(A_out,1)
    # On fait une pénalisation sur chaque nœud de Dirichlet
    for noeud in noeud_dirichlet
        # On fait une annulation de la ligne et de la colonne
        annuler_ligne_colonne(A_out,noeud)
        # On place 1 sur la diagonale pour imposer K[noeud,noeud]=1
        # Eviter que le système devienne singulier
        A_out[noeud,noeud]=1.0
        # Ensuite on met  b_out[noeud]=0 pour imposer u(noeud)=0
        b_out[noeud]=0.0
    end
    return A_out, b_out
end

function resoudre_systeme(A,b)
    N=size(A,1)  # Nombre d'inconnues qui est le nombre de ligne de la matrice
    # Si on a une petite matrice on utilise LU
    if N<5000
        return A\b  # Julia choisit LU ou QR selon la structure de A
    # Si on a une plus grande matrice symétrique définie positive, on utilise Cholesky
    else
        F=cholesky(A)
        return F\b
    end
end

function affichage_erreur(L, ux, NE_list)
    # Calcul des erreurs pour une source F1 et une solution exacte ux
    h, erreur = calcul_erreur(f1, ux, L, NE_list)
    
    # Tracé en échelle log-log
    plot(h, erreur, lw=2, xscale=:log10, yscale=:log10, xlabel="h (Taille du maillage)", ylabel="Erreur", title="Évolution de l'erreur en fonction de h (log-log)")
end

function calcul_erreur(f, ux, L, NE_list)
    # Calcul des pas de maillage h
    h = L ./ NE_list
    erreurs = Float64[]
    
    # Calcul de l'erreur pour chaque nombre d'éléments
    for NE in NE_list
        maillage_temporaire = "temp_maillage.txt"  # fichier temporaire de maillage
        
        # On génère le maillage temporaire
        maillage(L, NE, maillage_temporaire)
        
        # On lit le maillage temporaire
        liste_noeuds, liste_elements, noeud_dirichlet = lecture_maillage(maillage_temporaire)
        
        # Résolution du système pour obtenir la solution approximative
        solution_approchee = resoudre_systeme(A, b)  # Supposons que A et b sont calculés en fonction du maillage
        
        # Calcul de l'erreur (par exemple, erreur L2)
        erreur = calculer_erreur(ux, solution_approchee, liste_noeuds)
        push!(erreurs, erreur)
    end
    
    return h, erreurs
end

# Calcul de l'erreur (norme L2)
function calculer_erreur(ux, solution_approchee, liste_noeuds)
    erreur = 0.0
    # Suppose que ux est la fonction exacte et solution_approchee est la solution calculée
    for noeud in liste_noeuds
        # Calcul de l'erreur en un noeud
        erreur += abs(ux(noeud) - solution_approchee(noeud))^2  # Somme des erreurs au carré
    end
    return sqrt(erreur)
end


function maillage_2d(L,div,nom_fichier)
    out=open(nom_fichier,"w")
    NN=(div+1)^2
    NE=2*(div^2)
    # Nombre de noeuds NN et Nombre d'elements NE
    println(out,NN," ",NE)
    h=L/div
    noeud=1
    # Pour stocker chaque noeud dans conditions de Dirichlet
    noeud_cl_dirichlet=Int[]
    for j in 0:div
        for i in 0:div
            x=i*h
            y=j*h
            if x==0.0 || x==1.0 || y==0.0 || y==1.0
                push!(noeud_cl_dirichlet,noeud)
            end
            println(out,noeud," ",x," ", y)
            noeud+=1
        end
    end
    # on remarque que les coordonnées pour un carré sont:
    # en bas à gauche x,y : n1
    # en bas a droite x+1,y : n2
    # en haut à gauche x,y+1 : n3
    # en haut à droite x+1,y+1 : n4
    # donc on fait cela pour tout les carrées variant d'un certain nbr de div
    # triangle 1: n1 n2 n4
    # triangle 2: n2 n3 n4
    element = 1
    for y in 0:(div-1)
      for x in 0:(div-1)
        # Indices des 4 coins du carré local
        n1= y*(div+1)+x+1    # en bas à gauche
        n2=y*(div+1)+(x+1)+1    # en bas à droite
        n3=(y+1)*(div+1)+x+1   # en haut à gauche
        n4=(y+1)*(div+1)+(x+1)+1  # en haut à droit

        # Affichage du premier triangle
        println(out,element," ",n1, " ",n2, " ",n3)
        element += 1

        #Affichage du second triangle
        println(out,element," ",n2, " ",n4, " ",n3)
        element += 1
        end
    end
    # Affichage du nombre de noeud dans conditions de Dirichlet
    println(out,length(noeud_cl_dirichlet))
    # Affichage de tous ces noeuds
    for noeud in noeud_cl_dirichlet
        println(out,noeud)
    end
    close(out)
end


function imposition_cl_lim_alternative(A,b,liste_noeuds_cl_dirichlet)
    # On récupère le nombre de ligne de la matrice A
    n=size(A,1)
    # Créer vecteur pour stocker les indices à garder
    ind=[]
    # On récuypère les indices qui ne sont pas dans liste_noeuds_cl_dirichlet
    for i in 1:n
        if !(i in liste_noeuds_cl_dirichlet)
            push!(ind,i)
        end
    end
    # On créer la nouvelle matrice et le nouveau vecteur en utilisant avec seulemen les bons indices
    A_out=A[ind,ind]
    b_out=b[ind]
    return A_out,b_out
end

function reintegrer_uh(uh,liste_noeuds_cl_dirichlet,NN)
    # On initialise le vecteur total uh
    uh_tot=zeros(NN)
    j = 1
    # Remplit les valeurs manquantes en parcourant tous les noeuds
    for i in 1:NN
        if !(i in liste_noeuds_cl_dirichlet)
            uh_tot[i]=uh[j]
            j += 1
        end
    end
    return uh_tot
end



function affichage(A,b,liste_noeuds_cl_dirichlet)
    A_1,b_1=imposition_cl_lim_alternative(A,b,liste_noeuds_cl_dirichlet)
    uh=A_1\b_1
    uh_tot=reintegrer_uh(uh,liste_noeuds_cl_dirichlet,length(b_1))

    println("Matrice après suppression des lignes et des colonnes pour cl Dirichlet:")
    println(A_1)
    println("Vecteur du second membre après la modification:")
    println(b_1)
    println("La solution complète après la réintégration:")
    println(uh_tot)
end