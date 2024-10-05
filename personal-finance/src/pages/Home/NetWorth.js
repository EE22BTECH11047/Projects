import React, { Component } from 'react';

import Alert from '../../components/Alert';
import Storage from '../../components/Storage';
import Pie from '../../components/Pie';

import CONFIG from '../../config';

class NetWorth extends Component {
  constructor(props) {
    super(props);

    this.state = {
      all_items: Storage.getItems(null, props.items),
      assets: Storage.getItems(CONFIG.ITEM_TYPE.ASSET, props.items),
      liabilities: Storage.getItems(CONFIG.ITEM_TYPE.LIABILITY, props.items),
      income: Storage.getItems(CONFIG.ITEM_TYPE.INCOME, props.items),
      expenses: Storage.getItems(CONFIG.ITEM_TYPE.EXPENSE, props.items),
    };
  }

  render() {
    return (
      <React.Fragment>
        <div className="container">
          <br />

          <Alert items={this.state.all_items} />

          <br />

          <div className="row">
            <div className="col-sm">
              <Pie items={this.state.assets} title="Assets" />
            </div>

            <div className="col-sm"></div>

            <div className="col-sm">
              <Pie items={this.state.liabilities} title="Liabilities" />
            </div>
          </div>
        </div>
      </React.Fragment>
    );
  }
}

export default NetWorth;
