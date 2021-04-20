CHORaw = [50.14 6.16 43.22];

CHOPure = [CHORaw(1)/sum(CHORaw)...
           CHORaw(2)/sum(CHORaw)...
           CHORaw(3)/sum(CHORaw)]
   
sum(CHOPure)