

set.intersect <- function(set1,set2)
{
	set1 <- unique(set1);
	set2 <- unique(set2)
	set.cup <- c(set1,set2);
	set.cup.tab <- table(set.cup);
	set.cap <- as.numeric(names(set.cup.tab)[set.cup.tab==2]);
	return(set.cap);
}




set.compare <- function(set1,set2)
  {
    if(length(set1)!=length(set2)) return(0);
    if(length(set1)==length(set2))
      {
        return(sum(set1==set2)==length(set1));
      }
  }

