import React, { Component } from 'react';
import Highcharts from 'highcharts';
import highchartsExporting from 'highcharts/modules/exporting';
import HighchartsReact from 'highcharts-react-official';

highchartsExporting(Highcharts);

class Pie extends Component {
  constructor(props) {
    super(props);

    let data = [];

    if (props.items) {
      props.items.map(function (item) {
        data.push([item.name, item.amount]);

        return item;
      });
    }

    this.state = {
      options: {
        chart: {
          plotBackgroundColor: null,
          plotBorderWidth: null,
          plotShadow: false,
          backgroundColor: props.backgroundColor || 'transparent',
          type: 'pie',
        },
        title: {
          text: props.title,
        },
        tooltip: {
          pointFormat: '{series.name}: <b>{point.percentage:.1f}%</b>',
        },
        plotOptions: {
          pie: {
            allowPointSelect: true,
            cursor: 'pointer',
            dataLabels: {
              enabled: true,
              format: '<b>{point.name}</b>: {point.percentage:.1f} %',
            },
          },
        },
        series: [
          {
            name: 'Item',
            colorByPoint: true,
            data: data,
          },
        ],
      },
    };
  }

  render() {
    return (
      <HighchartsReact highcharts={Highcharts} options={this.state.options} />
    );
  }
}

export default Pie;
