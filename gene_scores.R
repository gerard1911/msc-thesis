# Create score matrix
gene.scores <- function(score){

scoremat = matrix(,nrow=5,ncol = 5)

scoremat[1,2] = score[1,1] + score[1,5] + score[1,6] + score[1,7] + score[1,11] + score[1,12] + score[1,13] + score[1,15]
scoremat[1,3] = score[1,2] + score[1,5] + score[1,8] + score[1,9] + score[1,11] + score[1,12] +  score[1,14] +score[1,15]
scoremat[1,4] = score[1,3] + score[1,6]  +score[1,8] + score[1,10] + score[1,11] + score[1,13] + score[1,14] + score[1,15]
scoremat[1,5] = score[1,4] + score[1,7]  +score[1,9] + score[1,10] + score[1,12] + score[1,13] +  score[1,14] + score[1,15]

scoremat[2,1] = score[2,1] + score[2,5] + score[2,6] + score[2,7] + score[2,11] + score[2,12] + score[2,13] + score[2,15]
scoremat[2,3] = score[2,2] + score[2,5] + score[2,8] + score[2,9] + score[2,11] + score[2,12] +  score[2,14] +score[2,15]
scoremat[2,4] = score[2,3] + score[2,6]  +score[2,8] + score[2,10] + score[2,11] + score[2,13] + score[2,14] + score[2,15]
scoremat[2,5] = score[2,4] + score[2,7]  +score[2,9] + score[2,10] + score[2,12] + score[2,13] +  score[2,14] + score[2,15]

scoremat[3,1] = score[3,1] + score[3,5] + score[3,6] + score[3,7] + score[3,11] + score[3,12] + score[3,13] + score[3,15]
scoremat[3,2] = score[3,2] + score[3,5] + score[3,8] + score[3,9] + score[3,11] + score[3,12] +  score[3,14] +score[3,15]
scoremat[3,4] = score[3,3] + score[3,6]  +score[3,8] + score[3,10] + score[3,11] + score[3,13] + score[3,14] + score[3,15]
scoremat[3,5] = score[3,4] + score[3,7]  +score[3,9] + score[3,10] + score[3,12] + score[3,13] +  score[3,14] + score[3,15]


scoremat[4,1] = score[4,1] + score[4,5] + score[4,6] + score[4,7] + score[4,11] + score[4,12] + score[4,13] + score[4,15]
scoremat[4,2] = score[4,2] + score[4,5] + score[4,8] + score[4,9] + score[4,11] + score[4,12] +  score[4,14] +score[4,15]
scoremat[4,3] = score[4,3] + score[4,6]  +score[4,8] + score[4,10] + score[4,11] + score[4,13] + score[4,14] + score[4,15]
scoremat[4,5] = score[4,4] + score[4,7]  +score[4,9] + score[4,10] + score[4,12] + score[4,13] +  score[4,14] + score[4,15]

scoremat[5,1] = score[5,1] + score[5,5] + score[5,6] + score[5,7] + score[5,11] + score[5,12] + score[5,13] + score[5,15]
scoremat[5,2] = score[5,2] + score[5,5] + score[5,8] + score[5,9] + score[5,11] + score[5,12] +  score[5,14] + score[5,15]
scoremat[5,3] = score[5,3] + score[5,6] + score[5,8] + score[5,10] + score[5,11] + score[5,13] + score[5,14] + score[5,15]
scoremat[5,4] = score[5,4] + score[5,7]  +score[5,9] + score[5,10] + score[5,12] + score[5,13] +  score[5,14] + score[5,15]
return(scoremat)
}