using DelimitedFiles
using LinearAlgebra
using SparseArrays
using Plots


############ETAPE 1##########################
function f1(x) # f qui a en entree x un nombre reel
    return 1.0 # sortie f(x) pre-remplit une valeur arbitraire
end

function f2(x) # f qui a en entree x un nombre reel
    return sin(x) # sortie f(x) pre-remplit une valeur arbitraire
end

############ETAPE 2##########################

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

#on utilise les meme parametres que mesh.txt
maillage(2,2,"maillage.txt")
# on voit bien que l'on a bien 11 noeuds , 10 elements
# l'affichage des noeuds est correcte on part de 0.0 et on va jusqu'a 1.0 en 11 neouds
# on affiche bien 10 elments ou chaque elements est représnté par 2 noeuds
# on a bien deux elemnts pour conditions de dirirchlet
# et les numeros sont biens 1 et 11


############ETAPE 3##########################

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

liste_noeuds,liste_elements,noeud_dirichlet=lecture_maillage("maillage.txt")
println("La liste des noeuds :")
println(liste_noeuds)
println("La liste des elements :")
println(liste_elements)
println("La liste des noeuds avec conditons de Dirichlet :")
println(noeud_dirichlet)
# ca affiche bien une liste avec la cordonnées de chacun des noeuds
# une liste avec la liste des noeuds de chaque elements
# et une liste avec les noeuds des conditions de dirichlet

############ETAPE 4##########################

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

# On fait le test sur l'element 1
Ke,Fe=tables_elementaires(1,liste_noeuds,liste_elements,f1)
println("Affichage de la matrice de rigidité Ke :")
println(Ke)
println("Affichage du second membre Fe :")
println(Fe)
# voir si les resultats correspondent bien a ce que on est censez avoir normalment oui


############ETAPE 5##########################

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

A_g,b_g=assembler_matrices(liste_noeuds,liste_elements,f1)
println("Affichage de la matrice globale :")
println(A_g)
println("Affichage du second membre globale :")
println(b_g)

############ETAPE 6##########################

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

A_cl_dirichlet,b_cl_dirichlet=imposer_cl_dirichlet(A_g,b_g,noeud_dirichlet)
println("Affichage de la matrice de rigidité globale avec cl de Dirichlet :")
println(A_cl_dirichlet)
println("Affichage du vecteur second membre avec cl de Dirichlet :")
println(b_cl_dirichlet)

############ETAPE 7##########################

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


###############ETAPE 8####################

# on définit la première solution exacte pour f1= 1.0
function ux_1(x, L)
    return .-(x.^2) ./ 2 .+ (L ./ 2) .* x
end

# on définit la deuxieme solution exacte pour f2=sin(x)
function ux_2(x, L)
    return sin.(x) .- (sin(L)/L) .* x
end

function affichage(f,L,NE,ux)
    # On génére le maillage
    maillage(L, NE, "maillage_1.txt")
    liste_noeuds, liste_elements, noeud_dirichlet = lecture_maillage("maillage_1.txt")
    # On fait l'assemblage
    A_g, b_g = assembler_matrices(liste_noeuds, liste_elements, f)
    # On appliquer les conditions de Dirichlet
    A_cl_dirichlet,b_cl_dirichlet = imposer_cl_dirichlet(A_g, b_g, noeud_dirichlet)
    # On résout le système
    u_h = resoudre_systeme(A_cl_dirichlet, b_cl_dirichlet)
    # On calcul la solution exacte sur un maillage pour pouvoir afficher
    x = range(0,L,200)
    f_res = ux(x,L)

    # On fait l'affichage final
    p=plot(x,f_res,label="Solution exacte",linewidth=2)
    plot!(liste_noeuds,u_h,label="Solution EF",linestyle=:dash, marker=:circle, linewidth=2)
    xlabel!("x")
    ylabel!("u(x)")
    title!("Comparaison entre la solution exacte et la solution par éléments finis")
    display(p)
end

# On créer les deux graphiques des elements finis
p1 = affichage(f1, 2.0, 2, ux_1)
p2 = affichage(f2, 10.0, 10, ux_2)

# On fait un affichage global
plot(p1, p2, layout=(1,2),size=(1500,400))

#################ETAPE 9#####################################################

function affichage_erreur(f, ux, L, NE_list)
    # Calcul des erreurs pour une source f et une solution exacte ux
    h, erreur = calcul_erreur(f, ux, L, NE_list)
    
    # Vérification de la validité des données avant de tracer
    println("h: ", h)
    println("erreur: ", erreur)
    
    # Tracé en échelle log-log
    plot(h, erreur, lw=2, xscale=:log10, yscale=:log10, xlabel="h (Taille du maillage)", ylabel="Erreur", title="Évolution de l'erreur en fonction de h (log-log)")
end

function calcul_erreur(f, ux, L, NE_list)
    # Calcul des pas de maillage h
    h = L ./ NE_list
    erreurs = Float64[]
    
    # Vérification si NE_list est vide
    if isempty(NE_list)
        println("La liste NE_list est vide.")
        return nothing, nothing  # Retourner Nothing pour éviter l'erreur
    end
    
    # Calcul de l'erreur pour chaque nombre d'éléments
    for NE in NE_list
        maillage_temporaire = "temp_maillage.txt"  # fichier temporaire de maillage
        
        # On génère le maillage temporaire
        maillage(L, NE, maillage_temporaire)
        
        # On lit le maillage temporaire
        liste_noeuds, liste_elements, noeud_dirichlet = lecture_maillage(maillage_temporaire)
        
        # Résolution du système pour obtenir la solution approximative
        A, b = assembler_matrice_et_vecteur(L, NE, f, liste_elements)
        solution_approchee = resoudre_systeme(A, b)  # Résolution du système
        
        # Calcul de l'erreur
        erreur = calculer_erreur(ux, solution_approchee, liste_noeuds)
        push!(erreurs, erreur)
    end
    
    # Vérification que des erreurs ont bien été calculées
    if isempty(erreurs)
        println("Aucune erreur calculée. Vérifiez les données.")
        return nothing, nothing  # Retourner Nothing si aucune erreur
    end
    
    return h, erreurs
end



L=10.0
NE_list=[5,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000]

# On créer les deux graphiques des erreurs
e1=affichage_erreur(f1,L,ux_1,NE_list)
e2=affichage_erreur(f1,L,ux_2,NE_list)

# On fait un affichage global
plot!(e1,e2,layout=(1,2),size=(1500,400))

#############ETAPE 10###############################

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

# réutilisation de A_g et b_g de l'étape 5
affichage(A_g,b_g,noeud_dirichlet)


 #############Étape 11##########################@
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

# en faisant une petite resolution de systeme on a que div = 10 dans exmeple donné mesh2d
div=10
L=1.0
maillage_2d(L,div,"maillage_2d.txt")

#########################Étape 12#####################################@
function lecture_maillage_2d(nom_fichier)
    #les positions des noeuds, les elements et les noeuds sur conditions de Dirichlet
    data=readdlm(nom_fichier) # le résultat est une matrice contenant toutes les valeurs du fichier
    # le nombre de noeuds
    NN=Int(data[1,1])
    # le nombre d'elements
    NE=Int(data[1,2])
    # les positions des noeuds
    liste_noeuds=hcat(data[2:NN+1,2],data[2:NN+1,3])
    # liste des noeuds de l'element
    liste_elements=data[NN+2:NN+1+NE,2:4]
    #Nombre Noeud sur condition de Dirichlet
    size_noeuds_dirichlet=data[NN+1+NE+1,1]
    # les noeuds sur conditions de Dirichlet
    noeud_dirichlet=data[NN+1+NE+2:NN+NE+2+size_noeuds_dirichlet,1]
    return liste_noeuds,liste_elements,noeud_dirichlet
end

liste_noeuds_2d,liste_elements_2d,noeud_dirichlet_2d=lecture_maillage_2d("mesh2D.msh")
println(liste_noeuds_2d)
println(liste_elements_2d)
println(noeud_dirichlet_2d)
# on a bien les bons resultats d'affiché



###################Étape 13############################@

# Utilisation de la fonction de lecture pour obtenir les données du maillage 2D
positions_noeud , liste_elements , dirichlet_noeuds = read_mesh2D("C:\\Users\\Adel\\Downloads\\mesh2D.msh")

# Fonction pour calculer la matrice de rigidité élémentaire d'un triangle P1
function table_elementaires2D(x1,y1,x2,y2,x3,y3,F)
aire = 0.5*abs.((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
b1,b2,b3 = y2-y3,y3-y1,y1-y2
c1,c2,c3 = x3-x2,x1-x3,x2-x1
# on utilise la forme de Ke retrouvé à la question 9 du td8
Ke = (1.0/(4.0 * aire))*[
b1*b1+c1*c1 b1*b2+c1*c2 b1*b3+c1*c3;
b2*b1+c2*c1 b2*b2+c2*c2 b2*b3+c2*c3;
b3*b1+c3*c1 b3*b2+c3*c2 b3*b3+c3*c3
]
F1,F2,F3 = F(x1) , F(x2) , F(x3)
Fe = (aire/3)*[F1,F2,F3] # on utilise la formule des moyennes
return Ke,Fe
end

function assemblage2D(positions_noeud, liste_elements, F)
# Calcul du nombre total de noeuds dans le maillage
NN = length(positions_noeud)
# Calcul du nombre d'éléments finis (chaque élément est défini par 3 noeuds)
M = Int(length(liste_elements) ÷ 3) # ÷ pour obtenir une division entière (nombre d'éléments)
# Initialisation de la matrice de rigidité globale (de taille NN x NN) et du vecteur de forces global (de taille NN)
K_global = zeros(NN, NN)
F_global = zeros(NN)
# Boucle sur chaque élément fini
for e in 1:M
# Récupération des indices des 3 noeuds de l'élément courant (e)
n1, n2, n3 = liste_elements[3*(e-1)+1], liste_elements[3*(e-1)+2], liste_elements[3*(e-1)+3]
# Récupération des coordonnées (x, y) des noeuds de l'élément courant
x1, y1 = positions_noeud[n1]
x2, y2 = positions_noeud[n2]
x3, y3 = positions_noeud[n3]
# Calcul de la matrice de rigidité élémentaire Ke et du vecteur de force élémentaire Fe
Ke, Fe = table_elementaires2D(x1, y1, x2, y2, x3, y3, F)
# Création d'une liste des noeuds de l'élément pour faciliter les boucles suivantes
noeuds = [n1, n2, n3]
# Boucle pour assembler les termes du vecteur de force global et de la matrice de rigidité globale
for i in 1:3
global_i = noeuds[i] # Récupère l'indice du nœud global i
# Assemblage du vecteur de force global
F_global[global_i] += Fe[i]
# Boucle interne pour assembler la matrice de rigidité globale
for j in 1:3
global_j = noeuds[j] # Récupère l'indice du nœud global j
# Assemblage des contributions dans la matrice de rigidité globale
K_global[global_i, global_j] += Ke[i, j]
end
end
end
# Retour de la matrice de rigidité globale et du vecteur de forces global assemblés
return K_global, F_global
end


# Etape 13 : visualisation de la solution obtenue

K_global , F_global = assemblage2D(positions_noeud,liste_elements,F1)
K_red , F_red = application_dirichlet2(K_global,F_global,dirichlet_noeuds)
uh_red = resolution_systeme(K_red,F_red)
uh = reintegration_dirichlet(uh_red,dirichlet_noeuds,NN)
X = [pos[1] for pos in positions_noeud]
Y = [pos[2] for pos in positions_noeud]
nx = length(unique(X))
ny = length(unique(Y))
Uh_matr = reshape(uh,nx,ny)'
plot(sort(unique(X)), sort(unique(Y)), Umat, st = :surface, xlabel="x", ylabel="y", zlabel="u(x,y)")
savefig("solution_2D_F1.png")


K_global , F_global = assemblage2D(positions_noeud,liste_elements,F2)
K_red , F_red = application_dirichlet2(K_global,F_global,dirichlet_noeuds)
uh_red = resolution_systeme(K_red,F_red)
uh = reintegration_dirichlet(uh_red,dirichlet_noeuds,NN)
X = [pos[1] for pos in positions_noeud]
Y = [pos[2] for pos in positions_noeud]
nx = length(unique(X))
ny = length(unique(Y))
Uh_matr = reshape(uh,nx,ny)'
plot(sort(unique(X)), sort(unique(Y)), Umat, st = :surface, xlabel="x", ylabel="y", zlabel="u(x,y)")
savefig("solution_2D_F2.png")


K_global , F_global = assemblage2D(positions_noeud,liste_elements,F3)
K_red , F_red = application_dirichlet2(K_global,F_global,dirichlet_noeuds)
uh_red = resolution_systeme(K_red,F_red)
uh = reintegration_dirichlet(uh_red,dirichlet_noeuds,NN)
X = [pos[1] for pos in positions_noeud]
Y = [pos[2] for pos in positions_noeud]
nx = length(unique(X))
ny = length(unique(Y))
Uh_matr = reshape(uh,nx,ny)'
plot(sort(unique(X)), sort(unique(Y)), Umat, st = :surface, xlabel="x", ylabel="y", zlabel="u(x,y)")
savefig("solution_2D_F3.png")