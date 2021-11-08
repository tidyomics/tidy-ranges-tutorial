library(dplyr)
library(tidyr)
library(randomNames)

m <- 17 # size of class
n <- 100 # size of population
nf <- 12 # number of friends

set.seed(5)
people <- randomNames(m + n,
                      name.order="first.last",
                      name.sep=" ")

pop <- data.frame(person=people[-(1:m)],
                  birthday=sample(31, n, TRUE))
friends <- data.frame(friend=sample(people, nf))
class_people <- people[1:m]

loop <- rep(1:m, each=nf)
idx <- as.vector(replicate(m, sample(n, nf)))

class <- data.frame(classmate=class_people[loop],
                    friend=pop$person[idx],
                    birthday=pop$birthday[idx])

class <- class %>%
  arrange(classmate) %>%
  mutate(classmate=factor(classmate))

shared_friends <- friends %>%
  inner_join(class, by="friend") %>%
  arrange(classmate, friend) %>%
  relocate(classmate)

shared_friends %>%
  group_by(classmate) %>%
  summarize(in_common=n(), median_bday=median(birthday)) %>%
  complete(classmate, fill=list(count=0))
