CHORaw = [49.59 6.06 44.08];

CHOPure = [CHORaw(1)/sum(CHORaw)...
           CHORaw(2)/sum(CHORaw)...
           CHORaw(3)/sum(CHORaw)]
   
sum(CHOPure)