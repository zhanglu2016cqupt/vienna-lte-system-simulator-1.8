classdef fifoQueue < handle
% Implements a FIFO queue with a circular buffer.
% (c) Josep Colom Ikuno, INTHFT, 2008

   properties
       % Where we store the data (numbers that will probably index
       % something)
       data
       % Where an insertion WILL be put
       insert_index
       % From where will we extract
       extract_index
       % Number of elements now in the queue
       size
       % circular buffer size
       buffer_size
   end

   methods
       % Class constructor. Specify a queue size_occupied big enough or you will
       % have 
       function obj = fifoQueue(queue_size)
           obj.data = zeros(1,queue_size);
           obj.insert_index  = queue_size;
           obj.extract_index = queue_size;
           obj.size = 0;
           obj.buffer_size = queue_size;
       end
       % Insert an element into the queue
       function insert(obj,data)
           % Check for overflow
           if obj.size>=obj.buffer_size
               error('Queue full');
           else
               % Insert data and move index
               obj.data(obj.insert_index) = data;
               obj.insert_index = obj.insert_index - 1;
               if obj.insert_index==0
                   obj.insert_index = obj.buffer_size;
               end
               obj.size = obj.size + 1;
           end
       end
       % Extract an element from the queue
       function data = extract(obj)
           % Check for overflow
           if obj.size<=0
               error('Queue empty');
           else
               % Insert data and move index
               data = obj.data(obj.extract_index);
               obj.extract_index = obj.extract_index - 1;
               if obj.extract_index==0
                   obj.extract_index = obj.buffer_size;
               end
               obj.size = obj.size - 1;
           end
       end
       % special function that will "delete" the entries that equal the
       % inputted value. Deleting means that it will put them to 0
       function delete(obj,data)
           obj.data = obj.data.*(obj.data~=data);
       end
       % Print some info
       function print(obj)
           fprintf('FIFO queue, size %d, %d occupied\n',obj.buffer_size,obj.size);
           if obj.size>0
               fprintf('Elements: ');
           end
           for i_=1:obj.size
               index = obj.extract_index-i_+1;
               if index<=0
                   index = index + obj.buffer_size;
               end
               fprintf('%d ',obj.data(index));
           end
           fprintf('\n');
       end
   end
end 
