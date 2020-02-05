// Molecular Orbital
( function() {

  var width = 150, height = 400;
  var svg = d3.select("#molecularorbital")
              .append("svg").attr("width",(width/800)*100+"%").attr("viewBox", [0,0, width, height]);
  var margin = {top: 20, right:40, bottom:30, left:40};

  //svg.style("background-color","pink");
  //svg.append("rect").attr("x", margin.right).attr("y", margin.top)
  //   .attr("width", width-(margin.left+margin.right))
  //   .attr("height", height-(margin.top+margin.bottom));

  var x = d3.scaleLinear()
            .domain([0, 1])
            .range([margin.left, width-margin.right])
            .nice();

  var y = d3.scaleLinear()
            .domain(d3.extent(data.mo_energy))
            .range([height-margin.bottom, margin.top])
            .nice();

  var vgrid = svg.append("g")
                 .selectAll("path")
                 .data(y.ticks())
                 .enter()
                 .append("line")
                 .attr("stroke-width", 1)
                 .attr("stroke-dasharray", "4")
                 .attr("stroke", "lightgray")
                 .attr("y1", d=>y(d))
                 .attr("y2", d=>y(d))
                 .attr("x1", margin.left)
                 .attr("x2", width-margin.right);

  irrepcolormap = {"A1":"red"};
  
  var mo = svg.append("g")
              .selectAll("path")
              .data(data.mo_energy)
              .enter().append("g").attr("class", "mo");
  mo.append("line")
    .attr("stroke-width", 1)
  //                 .attr("stroke-dasharray", "4")
    .attr("stroke", (d,i)=>(irrepcolormap[data.mo_irrep[i]] || "steelblue"))
    .attr("y1", d=>y(d))
    .attr("y2", d=>y(d))
    .attr("x1", margin.left)
    .attr("x2", width-margin.right);
  mo.append("text").attr("x", width-margin.right).attr("y", d=>y(d)).text((d,i)=>i);
  

  var yaxis = svg.append("g")
                 .call(d3.axisLeft(y))
                 .attr("transform", "translate(" + (margin.left) + ",0)" );

})();


// Hamiltonian
( function() {

  d3.select("#hamiltonian").append("p").text("" + (data.nterm.length-1) + " qubits (before tapering)");
  d3.select("#hamiltonian").append("p").text("" + (data.n_qubit) + " qubits (after tapering)");
  
  var width = 800, height = 400;
  var svg = d3.select("#hamiltonian")
              .append("svg").attr("viewBox", [0,0, width, height]);
  var margin = {top: 20, right:40, bottom:30, left:40};

  var x = d3.scaleLinear()
            .domain([-0.5, data.nterm.length-0.5])
            .range([margin.left, width-margin.right]);
            //.nice();

  var y = d3.scaleLinear()
            .domain(d3.extent(data.nterm))
            .range([height-margin.bottom, margin.top])
            .nice();

  var vgrid = svg.append("g")
                 .selectAll("path")
                 .data(y.ticks())
                 .enter()
                 .append("line")
                 .attr("stroke-width", 1)
                 .attr("stroke-dasharray", "4")
                 .attr("stroke", "lightgray")
                 .attr("y1", d=>y(d))
                 .attr("y2", d=>y(d))
                 .attr("x1", margin.left)
                 .attr("x2", width-margin.right);
  
  var hist_nterm = svg.append("g")
                      .selectAll("path")
                      .data(data.nterm)
                      .enter().append("g").attr("class", "mo");
   hist_nterm.append("rect")
    .attr("stroke-width", 1)
//    .attr("stroke", "steelblue")
    .attr("x", (d,i)=>x(i-0.4))
    .attr("y", d=>y(d))
    .attr("width", x(0.4)-x(0))
    .attr("height", d=>(y(0)-y(d)));

  var hist_nterm_tapered = svg.append("g")
                      .selectAll("path")
                      .data(data.nterm_tapered)
                      .enter().append("g").attr("class", "mo");
  
  hist_nterm_tapered.append("rect")
    .attr("stroke-width", 1)
//    .attr("stroke", "none")
    .attr("fill", "steelblue")
    .attr("x", (d,i)=>x(i))
    .attr("y", d=>y(d))
    .attr("width", x(0.4)-x(0))
    .attr("height", d=>(y(0)-y(d)));


  var xaxis = svg.append("g")
                 .call(d3.axisBottom(x))
                 .attr("transform", "translate(0," + (height-margin.bottom) + ")" );
  var yaxis = svg.append("g")
                 .call(d3.axisLeft(y))
                 .attr("transform", "translate(" + (margin.left) + ",0)" );

})();


// Simulation
( function() {

  var width = 800, height = 400;
  var svg = d3.select("#simulation")
              .append("svg").attr("viewBox", [0,0, width, height]);
  var margin = {top: 20, right:40, bottom:30, left:40};

  var x = d3.scaleLinear()
            .domain([-0.5, data.ngate.length-0.5])
            .range([margin.left, width-margin.right]);
            //.nice();

  var y = d3.scaleLinear()
            .domain(d3.extent(data.ngate))
            .range([height-margin.bottom, margin.top])
            .nice();

  var vgrid = svg.append("g")
                 .selectAll("path")
                 .data(y.ticks())
                 .enter()
                 .append("line")
                 .attr("stroke-width", 1)
                 .attr("stroke-dasharray", "4")
                 .attr("stroke", "lightgray")
                 .attr("y1", d=>y(d))
                 .attr("y2", d=>y(d))
                 .attr("x1", margin.left)
                 .attr("x2", width-margin.right);
  
  var hist_ngate = svg.append("g")
                      .selectAll("path")
                      .data(data.ngate)
                      .enter().append("g").attr("class", "mo");
   hist_ngate.append("rect")
    .attr("stroke-width", 1)
  //    .attr("stroke", "steelblue")
    .attr("fill", "steelblue")  
    .attr("x", (d,i)=>x(i-0.4))
    .attr("y", d=>y(d))
    .attr("width", x(0.8)-x(0))
    .attr("height", d=>(y(0)-y(d)));

  var xaxis = svg.append("g")
                 .call(d3.axisBottom(x))
                 .attr("transform", "translate(0," + (height-margin.bottom) + ")" );
  var yaxis = svg.append("g")
                 .call(d3.axisLeft(y))
                 .attr("transform", "translate(" + (margin.left) + ",0)" );

  var time = d3.select("#simulation").append("div").style("font-size", "150%");
  time.append("p").attr("viewBox", [0,0, width, height]).text("Total Simulation Time : " + data.time_sim + " second");
  time.append("p").attr("viewBox", [0,0, width, height]).text("Num of trotter steps : " + data.n_trott_step);
  time.append("p").attr("viewBox", [0,0, width, height]).text("Simulation time per step : " + data.time_sim/data.n_trott_step + " second");

})();

// Correlation Function
( function() {

  //d3.select("body").style("margin", 0);
  
  var width = 800, height = 400;
  var svg = d3.select("#correlationfunction")
              //.style("background-color", "rgb(240, 255, 255)").style("padding-left", 4)
              .append("svg").attr("viewBox", [0,0, width, height]);
  var margin = {top: 20, right:40, bottom:30, left:40};

  /*
  svg.append("rect")
     .attr("x", 0)
     .attr("y", 0)
     .attr("width", width)
     .attr("height", height)
     .attr("fill","white");

  svg.append("rect")
     .attr("x", margin.left)
     .attr("y", margin.top)
     .attr("width", width-(margin.right + margin.left))
     .attr("height", height-(margin.top+margin.bottom))
     .attr("fill","white");
  */

  var time = data.corr_re_exact.map((d,i)=>i*data.dt);
  console.log(time);

  var x = d3.scaleLinear()
            .domain(d3.extent(time))
            .range([margin.left, width-margin.right])
            .nice();

  var y = d3.scaleLinear()
            .domain([-1,1])
            .range([height-margin.bottom, margin.top])
            .nice();

  var vgrid = svg.append("g")
                 .selectAll("path")
                 .data(x.ticks())
                 .enter()
                 .append("line")
                 .attr("stroke-width", 1)
                 .attr("stroke-dasharray", "4")
                 .attr("stroke", "lightgray")
                 .attr("x1", d=>x(d))
                 .attr("x2", d=>x(d))
                 .attr("y1", height-margin.bottom)
                 .attr("y2", margin.top);

  var line = d3.line().x((d,i)=>x(time[i])).y(d=>y(d));

  svg.append("path")
     .datum(data.corr_re_trott)
     .attr("fill", "none")
     .attr("stroke", "black")
     .attr("stroke-width", 1.5)
     .attr("d", line);
  
  svg.append("path")
     .datum(data.corr_im_trott)
     .attr("fill", "none")
     .attr("stroke", "black")
     .attr("stroke-width", 1.5)
     .attr("stroke-dasharray", '4')
     .attr("d", line);

  svg.append("path")
     .datum(data.corr_re_exact)
     .attr("fill", "none")
     .attr("stroke", "steelblue")
     .attr("stroke-width", 1.5)
     .attr("d", line);

  svg.append("path")
     .datum(data.corr_im_exact)
     .attr("fill", "none")
     .attr("stroke", "steelblue")
     .attr("stroke-width", 1.5)
     .attr("stroke-dasharray", '4')
     .attr("d", line);

  
  var xaxis = svg.append("g")
                 .call(d3.axisBottom(x))
                 .attr("transform", "translate(0," + y(0) + ")" );
                 //.attr("transform", "translate(0," + (height-margin.bottom) + ")" );

  var yaxis = svg.append("g")
                 .call(d3.axisLeft(y))
                 .attr("transform", "translate(" + (margin.left) + ",0)" );

})();

// Spectrum
( function() {

  var width = 800, height = 400;
  var svg = d3.select("#spectrum").append("svg").attr("viewBox", [0,0, width, height]);
  var margin = {top: 20, right:40, bottom:30, left:40};

  
/*
  svg.append("rect")
     .attr("x", 0)
     .attr("y", 0)
     .attr("width", width)
     .attr("height", height)
     .attr("fill","white");

  svg.append("rect")
     .attr("x", margin.left)
     .attr("y", margin.top)
     .attr("width", width-(margin.right + margin.left))
     .attr("height", height-(margin.top+margin.bottom))
     .attr("fill","white");
  */

  var x = d3.scaleLinear()
            .domain(d3.extent(data.energy_range))
            .range([margin.left, width-margin.right]);
            //.nice();

  var y = d3.scaleLinear()
            .domain([0, d3.max(data.spectrum_exact)*1.02])
            .range([height-margin.bottom, margin.top]);
            //.nice();

  var y_prony = d3.scaleLinear()
                  .domain([0,1])
                  .range([height-margin.bottom, margin.top])
                  .nice();

  var vgrid = svg.append("g")
                 .selectAll("path")
                 .data(x.ticks())
                 .enter()
                 .append("line")
                 .attr("stroke-width", 1)
                 .attr("stroke-dasharray", "4")
                 .attr("stroke", "lightgray")
                 .attr("x1", d=>x(d))
                 .attr("x2", d=>x(d))
                 .attr("y1", height-margin.bottom)
                 .attr("y2", margin.top);

  var eigvals = svg.append("g")
                   .selectAll("path")
                   .data(data.energy)
                   .enter()
                   .append("line")
                   .attr("stroke-width", 1.5)
                   .attr("stroke", "steelblue")
                   .attr("x1", d=>x(d))
                   .attr("x2", d=>x(d))
                   .attr("y1", height-margin.bottom)
                   .attr("y2", margin.top);

  /*
  var prony_exact = svg.append("g").selectAll("xxx").data(data.prony_p_exact).enter().append("circle")
                       .attr("cx", d=>x((-d-2*Math.PI)/data.dt) )
                       .attr("cy", (d,i)=>y_prony(data.prony_A_exact[i]))
                       .attr("r", 3)
                       .attr("stroke","red")
                       .attr("stroke-width",1)
                       .attr("fill","none");

  var prony_trott = svg.append("g").selectAll("xxx").data(data.prony_p_trott).enter().append("circle")
                       .attr("cx", d=>x((-d-2*Math.PI)/data.dt) )
                       .attr("cy", (d,i)=>y_prony(data.prony_A_trott[i]))
                       .attr("r", 3)
                       .attr("stroke","green")
                       .attr("stroke-width",1)
                       .attr("fill","none");
  */

  var prony_exact = svg.append("g").selectAll("xxx").data(data.prony_e_exact).enter().append("circle")
                       .attr("cx", d=>x(d))
                       .attr("cy", (d,i)=>y_prony(data.prony_a_exact[i]))
                       .attr("r", 3)
                       .attr("stroke","black")
                       .attr("stroke-width",1)
                       .attr("fill", "none");

  var prony_trott = svg.append("g").selectAll("xxx").data(data.prony_e_trott).enter().append("circle")
                       .attr("cx", d=>x(d))
                       .attr("cy", (d,i)=>y_prony(data.prony_a_trott[i]))
                       .attr("r", 3)
                       .attr("stroke","blue")
                       .attr("stroke-width",1)
                       .attr("fill", "none");


  var freq_vs_amp_exact = [];
  for(var i=0; i<data.energy_range.length; ++i) {
    freq_vs_amp_exact.push([data.energy_range[i], data.spectrum_exact[i]]);
  }

  var freq_vs_amp_trott = [];
  for(var i=0; i<data.energy_range.length; ++i) {
    freq_vs_amp_trott.push([data.energy_range[i], data.spectrum_trott[i]]);
  }
  
  var line = d3.line().x(d=>x(d[0])).y(d=>y(d[1]));

  svg.append("path")
     .datum(freq_vs_amp_exact)
     .attr("fill", "none")
     .attr("stroke", "black")
     .attr("stroke-width", 1.5)
     .attr("d", line);

  svg.append("path")
     .datum(freq_vs_amp_trott)
     .attr("fill", "none")
     .attr("stroke", "blue")
     .attr("stroke-width", 1.5)
     .attr("d", line);
  
  var xaxis = svg.append("g")
                 .call(d3.axisBottom(x))
                 .attr("transform", "translate(0," + (height-margin.bottom) + ")" );

  var yaxis = svg.append("g")
                 .call(d3.axisLeft(y))
                 .attr("transform", "translate(" + (margin.left) + ",0)" );

  var yaxis = svg.append("g")
                 .call(d3.axisRight(y_prony))
                 .attr("transform", "translate(" + (width-margin.right) + ",0)" );
})();
