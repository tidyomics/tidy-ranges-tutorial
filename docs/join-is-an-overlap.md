# Join is an overlap

Objective: learn that a "join" is like an overlap.

We start with a quick example explaining why we use functions called
`join`.


```r
library(dplyr)
library(tidyr)
library(randomNames)
```

Let's set up a simulation where we have 17 classmates (not counting
ourselves) in a school of 118. Suppose every member of the class has
12 friends outside of class. 

We want to find out, for each classmate, how many friends we have in
common and also, of our shared friends, what is a typical birthday.
Let's define typical as the median birthday of our shared friends.


```r
m <- 17 # size of class
n <- 100 # size of others in school
nf <- 12 # number of friends outside class
set.seed(5)
people <- randomNames(m + n,
                      name.order="first.last",
                      name.sep=" ")
```

We define the population of potential friends (those outside the
class) as `pop`, and our 12 friends are in `friends`. Finally, we
define the people that are in our class as `class_people`.


```r
pop <- data.frame(person=people[-(1:m)],
                  birthday=sample(31, n, TRUE))
friends <- data.frame(friend=sample(people, nf))
class_people <- people[1:m]
```

The following sets up a data.frame, where each row gives, for a given
classmate, their friends, and the friends' birthdays.


```r
loop <- rep(1:m, each=nf)
idx <- as.vector(replicate(m, sample(n, nf)))
class <- data.frame(classmate=class_people[loop],
                    friend=pop$person[idx],
                    birthday=pop$birthday[idx])
```

Sort this by classmate alphabetically, and make classmate into a
factor. This last step is important, as it will help us to keep track
of the classmates for whom we share _no_ friends outside of class.


```r
class <- class %>%
  arrange(classmate) %>%
  mutate(classmate=factor(classmate))
```

We perform an `inner_join` by `"friend"`. This just means we look for
classmates where we have overlapping friends, and we drop the rows
where we don't share any friends. "Inner" refers to the fact that we
are keeping the overlap in the _inside_ of two intersecting circles.
Note that the join operation brings along the metadata (extra data)
about the friends' birthdays.


```r
shared_friends <- friends %>%
  inner_join(class, by="friend") %>%
  arrange(classmate, friend) %>%
  relocate(classmate) # classmate to 1st column
shared_friends
```

```
##                classmate             friend birthday
## 1          Brandon Jones Christopher Rivera       13
## 2          Brandon Jones       Jose Jimenez       26
## 3          Brandon Jones    Maria Hernandez       15
## 4          Brandon Jones         Tuli Hoang        7
## 5     Christopher Herron         Tuli Hoang        7
## 6     Christopher Herron  William Steinbach       15
## 7         Haley Polhamus   Cellene Millhone       23
## 8         Haley Polhamus  William Steinbach       15
## 9  Juan Villegas Cabrera     Alyssa Kinanee       30
## 10 Juan Villegas Cabrera       Malik Gammon        5
## 11 Juan Villegas Cabrera     Tajhae Bohanna       31
## 12     Maisara el-Arshad     Kaylyn Judkins        2
## 13     Maisara el-Arshad    Maria Hernandez       15
## 14     Maisara el-Arshad     Tajhae Bohanna       31
## 15        Michael Mcgill Christopher Rivera       13
## 16        Michael Mcgill    Maria Hernandez       15
## 17        Noah Pettinger     Kaylyn Judkins        2
## 18        Noah Pettinger     Tajhae Bohanna       31
## 19         Orion Villani       Malik Gammon        5
## 20         Shane Ranaldi     Alyssa Kinanee       30
## 21         Shane Ranaldi  William Steinbach       15
## 22   Stephanie Hernandez Asmaa el-Abdelnour       29
## 23   Zachary Roe-Huffman    Maria Hernandez       15
```

Lastly, we perform some summarization: compute the number of friends
in common with `n()` and the median birthday of shared friends. The
`complete` call at the end fills in 0 for those classmates for
whom we share no friends (here, the use of `factor` earlier becomes
relevant). We can choose which columns to fill in, and what value to 
add.


```r
shared_friends %>%
  group_by(classmate) %>%
  summarize(in_common=n(), median_bday=median(birthday)) %>%
  complete(classmate, fill=list(in_common=0,median_bday=-1))
```

```
## # A tibble: 17 × 3
##    classmate             in_common median_bday
##    <fct>                     <int>       <dbl>
##  1 Aaliyah Minter                0        -1  
##  2 Brandon Jones                 4        14  
##  3 Christopher Herron            2        11  
##  4 Collin Leon                   0        -1  
##  5 Haley Polhamus                2        19  
##  6 Juan Villegas Cabrera         3        30  
##  7 Khaalid el-Ammar              0        -1  
##  8 Kianna Mcalevy                0        -1  
##  9 Maazin al-Ismael              0        -1  
## 10 Maisara el-Arshad             3        15  
## 11 Michael Mcgill                2        14  
## 12 Noah Pettinger                2        16.5
## 13 Orion Villani                 1         5  
## 14 Shane Ranaldi                 2        22.5
## 15 Stephanie Hernandez           1        29  
## 16 Xavier Urueta                 0        -1  
## 17 Zachary Roe-Huffman           1        15
```

