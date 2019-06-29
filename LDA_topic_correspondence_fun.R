# Computing a one-to-one topic correspondence between j_real and j_select: 
# every topic in j_select is associated with a single topic in j_real,
# the best matching among the topics that are not assigned yet. 

LDA_topic_correspondence_fun <- function(KL_topic_comparison,maxmatching,greedymatching){
  
  nb_topics = nrow(KL_topic_comparison)
  Topic_correspondence = vector(length=nb_topics,mode="numeric")
  
  if (maxmatching)
  {
    KL_topic_comparison.coor = t(as.matrix(expand.grid(1:nb_topics,1:nb_topics)))
    # Each pair in KL_topic_comparison.coor corresponds to a pair (k,k1), 
    # the k1 elements are relabelled between nb_topics+1 and 2*nb_topics in topic_graph
    topic_graph = KL_topic_comparison.coor
    topic_graph[2,] = topic_graph[2,] + nb_topics
    topic_graph = igraph::graph(as.vector(topic_graph), directed = F)
    # The edges are weighted using the opposite of the KL_topic_comparison matrix, because the algorithm maximizes the weights instead of minimizing them
    # Weights are also transformed into integers for handling by the underlying C code in max_bipartite_match
    topic_graph_weights = vector(length = nb_topics*nb_topics, mode= "numeric")
    KL_topic_comparison = round((max(KL_topic_comparison) - KL_topic_comparison)*100000)
    # Replacing the 0 value by 1, so as to keep all topics connected (otherwise the link with weight 0 cannot be assigned as a best match by the algo)
    KL_topic_comparison[which(KL_topic_comparison == 0)] = 1
    for (i in 1:(nb_topics*nb_topics))
      topic_graph_weights[i] = KL_topic_comparison[KL_topic_comparison.coor[1,i],KL_topic_comparison.coor[2,i]]
    igraph::E(topic_graph)$weight = topic_graph_weights
    # Matching is done using the Hungarian algorithm by calling the maxmatching function of package "maxmatching",
    # which merely calls the function max_bipartite_match of package "igraph"
    matching = maxmatching(topic_graph, weighted = T) 
    for (k in 1:nb_topics)
      Topic_correspondence[k] = matching$matching[k] - nb_topics 
    # Topic_correspondence1 = vector(length=nb_topics,mode="numeric")
    # for (k in (nb_topics+1):(2*nb_topics))
    #   Topic_correspondence1[matching$matching[k]] = k - nb_topics 
  } else if (greedymatching)
  {
    i = 1
    nb_topics_assigned = 0
    k_vect = vector(length=nb_topics,mode="numeric")
    k1_vect = vector(length=nb_topics,mode="numeric")
    sorted_KL_topic_comparison = sort.int(KL_topic_comparison,index.return=T)
    while (nb_topics_assigned < nb_topics)
    {
      k = sorted_KL_topic_comparison$ix[i] %% nb_topics 
      k1 = sorted_KL_topic_comparison$ix[i] %/% nb_topics
      if (k==0)
        k = nb_topics
      else 
        k1 = k1 + 1
      if (!any(k_vect == k) && !any(k1_vect == k1))
      {
        Topic_correspondence[k] = k1
        nb_topics_assigned = nb_topics_assigned+1
        k_vect[nb_topics_assigned] = k
        k1_vect[nb_topics_assigned] = k1
      }
      i = i+1
    }
  }
  
  return(Topic_correspondence)
}