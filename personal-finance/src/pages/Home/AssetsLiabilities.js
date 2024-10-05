import React, { Component } from 'react';

import Assets from '../../components/Assets';
import Liabilities from '../../components/Liabilities';
import Storage from '../../components/Storage';

import CONFIG from '../../config';

class IncomeExpenses extends Component {
  constructor(props) {
    super(props);

    this.state = {
      all_items: Storage.getItems(null, props.items),
      assets: Storage.getItems(CONFIG.ITEM_TYPE.ASSET, props.items),
      liabilities: Storage.getItems(CONFIG.ITEM_TYPE.LIABILITY, props.items),
    };
  }

  render() {
    return (
      <React.Fragment>
        <div className="container">
          <br />

          <div className="row">
            <div className="col-sm">
              <Assets items={this.state.assets} />
            </div>

            <div className="col-sm">
              <Liabilities items={this.state.liabilities} />
            </div>
          </div>
        </div>
      </React.Fragment>
    );
  }
}

export default IncomeExpenses;
