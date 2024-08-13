#!/bin/R Scripts 
# Date: 8/8/2024
# Author: Kasonde Chewe 
# > Description: Scripts below are helps in the data extraction, processing and 
# >              visualization. See function definitions for details


# Function to install and load packages
install_and_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, repos = "https://cloud.r-project.org/")
    library(package, character.only = TRUE)
  }
}


# Function to access an XML node and its subnode by name
get_node_value_by_name <- function(root, node_name, subnode_name) {
  # XPath to get the specific node and subnode
  xpath <- paste0("//", node_name, "/", subnode_name)
  
  # Get the node set using XPath
  nodes <- getNodeSet(root, xpath)
  
  # Return the value of the first matched node (if any)
  if (length(nodes) > 0) {
    return(xmlValue(nodes[[1]]))
  } else {
    return(NULL)
  }
}


# Function to extract values from a node
extract_values <- function(node, fields) {
  values <- sapply(fields, function(field) {
    value_node <- getNodeSet(node, paste0(".//", field))
    if (length(value_node) > 0) {
      return(xmlValue(value_node[[1]]))
    } else {
      return(NA)
    }
  })
  return(values)
}

# Function to recursively list all nodes and their attributes
list_nodes_to_file <- function(node, file_conn, indent = "") {
  # Write the current node name to file
  writeLines(paste0(indent, xmlName(node)), file_conn)
  
  # If the node has attributes, write them to file
  if (length(xmlAttrs(node)) > 0) {
    writeLines(paste0(indent, "  Attributes:"), file_conn)
    for (attr_name in names(xmlAttrs(node))) {
      writeLines(paste0(indent, "    -", attr_name, ":", xmlAttrs(node)[[attr_name]]), file_conn)
    }
  }
  
  # If the node has children, recursively list them to file
  if (length(xmlChildren(node)) > 0) {
    for (child in xmlChildren(node)) {
      if (inherits(child, "XMLInternalElementNode")) {
        list_nodes_to_file(child, file_conn, paste0(indent, "  "))
      }
    }
  }
}