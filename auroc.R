auroc <- function(score.matrix,relation.matrix){

EMscore_v = as.vector(score.matrix) 
relation_v = as.vector(relation.matrix)
relation_v = relation_v[!is.na(EMscore_v)]
EMscore_v = EMscore_v[!is.na(EMscore_v)]
order_score = order(EMscore_v)
relation_v = relation_v[order(EMscore_v, partial = relation_v, decreasing = TRUE)]
EMscore_v = sort(EMscore_v, decreasing = TRUE)

# Create ROC vector
tp = c()
fp = c()

tp[1] = 0
fp[1] = 0
k1 = 0
k2 = 0
auroc_score = 0
i=1
j=1
while(i <=length(EMscore_v)){
	if(relation_v[i]==TRUE && sum(EMscore_v==EMscore_v[i]) == 1 ){
		tp[j+1] = tp[j] + 1/ sum(relation_v)
		fp[j+1] = fp[j]	
		}else if(relation_v[i]==FALSE && sum(EMscore_v==EMscore_v[i]) == 1 ){
		tp[j+1] = tp[j] 
		fp[j+1] = fp[j] + 1 / sum(relation_v==FALSE)
		k1=k2
		k2 = fp[j+1]
		auroc_score = auroc_score + tp[j+1]*(k2-k1)
		}else if(relation_v[i]==TRUE && sum(EMscore_v==EMscore_v[i]) != 1 ){
			tp[j+1] = tp[j] + sum(relation_v[EMscore_v==EMscore_v[i]]==TRUE)/sum(relation_v)
			fp[j+1] = fp[j] + sum(relation_v[EMscore_v==EMscore_v[i]] == FALSE)/sum(relation_v==FALSE)
			k1=k2
			k2 = fp[j+1]
			auroc_score = auroc_score + 0.5*(tp[j+1]-tp[j])*(k2-k1) + tp[j]*(k2-k1) 

		}else if(relation_v[i]==FALSE && sum(EMscore_v==EMscore_v[i]) != 1 ){
			tp[j+1] = tp[j] + sum(relation_v[EMscore_v==EMscore_v[i]]==TRUE) /sum(relation_v)
			fp[j+1] = fp[j] + sum(relation_v[EMscore_v==EMscore_v[i]] == FALSE)/sum(relation_v==FALSE)
			k1=k2
			k2 = fp[j+1]
			auroc_score = auroc_score + 0.5*(tp[j+1]-tp[j])*(k2-k1) + tp[j]*(k2-k1) 
		}
		i = i + sum(EMscore_v==EMscore_v[i])
		j = j+1
}
plot(fp,tp, main = 'ROC-curve', xlab = 'false positive', ylab = 'true positive')
lines(fp,tp)
return(list(AUROC = auroc_score, ROC =  cbind(fp,tp)))
}