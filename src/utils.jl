# Function to get all subfolders
function get_subfolder(folder)
    # Get all subfolders
    paths = string.(folder,"/",readdir(folder))
    what = isdir.(paths) .& .!(occursin.("transport_properties",paths)) .& .!(startswith.(readdir(folder),"X"))
    subf = paths[what]
end
